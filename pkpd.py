"""药动学（PK）与量-效-时（PD）模型：模拟与拟合。

模型
----
量-效（Hill / Emax 方程）：
    E(D) = E0 + Emax * D^h / (EC50^h + D^h)

一室一级吸收药动学（血管外给药）：
    C(t) = (F·D·ka) / (V·(ka - ke)) · (e^(-ke·t) - e^(-ka·t))
    ka ≈ ke 时退化为 (F·D·ka / V) · t · e^(-ka·t)

时-效（PK 驱动）：把 C(t) 代入 Hill 方程，
    E(t) = E0 + Emax · C(t)^h / (EC50^h + C(t)^h)

拟合用 scipy.optimize.curve_fit，输入为实验测得的 (剂量, 效应) 或 (时间, 浓度)。
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


# ----------------------------------------------------------------- 模型公式
def hill_emax(dose, e0: float, emax: float, ec50: float, h: float):
    """量-效 Hill 方程。dose 为标量或数组。"""
    dose = np.asarray(dose, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        value = e0 + emax * np.power(dose, h) / (np.power(ec50, h) + np.power(dose, h))
    return np.where(dose == 0, e0, value)


def pk_concentration(t, dose: float, f: float, v: float,
                     ka: float, ke: float):
    """一室一级吸收的血药浓度 C(t)。"""
    t = np.asarray(t, dtype=float)
    if abs(ka - ke) < 1e-6:  # 退化情形，避免除零
        return (f * dose * ka / v) * t * np.exp(-ka * t)
    return (f * dose * ka) / (v * (ka - ke)) * (np.exp(-ke * t) - np.exp(-ka * t))


def effect_time(t, dose: float, f: float, v: float, ka: float, ke: float,
                e0: float, emax: float, ec50: float, h: float):
    """PK 驱动的时-效曲线。"""
    return hill_emax(pk_concentration(t, dose, f, v, ka, ke), e0, emax, ec50, h)


# ----------------------------------------------------------------- 曲线模拟
def simulate_pk(dose: float, f: float, v: float, ka: float, ke: float,
                t_max: float = 24.0, n: int = 241) -> pd.DataFrame:
    t = np.linspace(0.0, t_max, n)
    c = pk_concentration(t, dose, f, v, ka, ke)
    cmax, tmax = float(c.max()), float(t[np.argmax(c)])
    return pd.DataFrame({"时间 (h)": t, "血药浓度 (mg/L)": c}), cmax, tmax


def simulate_dose_response(e0: float, emax: float, ec50: float, h: float,
                           d_max: float = 200.0, n: int = 200) -> pd.DataFrame:
    dose = np.linspace(0.0, d_max, n)
    return pd.DataFrame({"剂量 (mg/kg)": dose,
                         "效应 (%)": hill_emax(dose, e0, emax, ec50, h)})


def simulate_effect_time(dose: float, f: float, v: float, ka: float, ke: float,
                         e0: float, emax: float, ec50: float, h: float,
                         t_max: float = 24.0, n: int = 241) -> pd.DataFrame:
    t = np.linspace(0.0, t_max, n)
    return pd.DataFrame({"时间 (h)": t, "效应 (%)": effect_time(
        t, dose, f, v, ka, ke, e0, emax, ec50, h)})


# ----------------------------------------------------------------- 实验拟合
def _r_squared(y_true, y_pred) -> float:
    y_true = np.asarray(y_true, dtype=float)
    ss_res = float(np.sum((y_true - y_pred) ** 2))
    ss_tot = float(np.sum((y_true - np.mean(y_true)) ** 2))
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")


def fit_dose_response(dose, effect) -> Tuple[Optional[Dict[str, float]], float]:
    """拟合 (剂量, 效应) 数据，返回 (参数字典, R²)。"""
    dose = np.asarray(dose, dtype=float)
    effect = np.asarray(effect, dtype=float)
    if len(dose) < 4:
        return None, float("nan")
    try:
        popt, _ = curve_fit(
            hill_emax, dose, effect,
            p0=[min(effect), max(effect) - min(effect), np.median(dose[dose > 0]), 1.0],
            bounds=([-np.inf, 0, 1e-6, 0.1], [np.inf, np.inf, np.inf, 10]),
            maxfev=20000,
        )
    except Exception:
        return None, float("nan")
    params = {"E0": popt[0], "Emax": popt[1], "EC50": popt[2], "h": popt[3]}
    return params, _r_squared(effect, hill_emax(dose, *popt))


def fit_pk(t, conc) -> Tuple[Optional[Dict[str, float]], float]:
    """拟合 (时间, 浓度) 数据，返回 (参数字典, R²)；A = F·D·ka/(V(ka-ke))。"""
    t = np.asarray(t, dtype=float)
    conc = np.asarray(conc, dtype=float)
    if len(t) < 4:
        return None, float("nan")

    def model(tt, a, ka, ke):
        if abs(ka - ke) < 1e-6:
            return a * ka * tt * np.exp(-ka * tt)
        return a * (np.exp(-ke * tt) - np.exp(-ka * tt))

    try:
        popt, _ = curve_fit(model, t, conc, p0=[max(conc), 1.0, 0.1],
                            bounds=([0, 1e-4, 1e-4], [np.inf, 20, 20]), maxfev=20000)
    except Exception:
        return None, float("nan")
    params = {"A": popt[0], "ka": popt[1], "ke": popt[2],
              "t1/2(吸收)": np.log(2) / popt[1], "t1/2(消除)": np.log(2) / popt[2]}
    return params, _r_squared(conc, model(t, *popt))


def pk_metrics(ka: float, ke: float, dose: float, f: float, v: float) -> Dict[str, float]:
    """常用药动学指标：达峰时间、半衰期、表观分布容积换算的 AUC（梯形近似）。"""
    frame, cmax, tmax = simulate_pk(dose, f, v, ka, ke)
    auc = float(np.trapezoid(frame["血药浓度 (mg/L)"], frame["时间 (h)"]))
    return {"Cmax": cmax, "Tmax": tmax, "AUC0-t": auc,
            "t1/2(ka)": float(np.log(2) / ka), "t1/2(ke)": float(np.log(2) / ke)}
