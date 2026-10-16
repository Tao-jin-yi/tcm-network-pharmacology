"""把外部的中药-成分数据表导入本项目（用于录入 TCMSP 全部 499 味中药）。

支持 CSV / Excel。表头可以是中文或英文，脚本会自动识别：

    中药名         herb / 中药 / 药材 / drug
    成分编号       mol_id / MOLID / 成分ID / molecule_id
    成分名         name / 成分 / 成分名 / molecule_name
    英文名         en / 英文名 / english_name
    口服生物利用度  ob / OB / 口服生物利用度 / bioavailability
    类药性         dl / DL / 类药性 / drug_likeness
    靶点           targets / 靶点 / gene / target   （多个用 ; 或 , 分隔）

用法：
    python import_tcmsp.py TCMSP导出.csv
    python import_tcmsp.py TCMSP导出.xlsx --replace   # 丢弃内置示例数据，只用导入的数据

脚本会生成 tcmsp_dataset.json；随后运行 python data_prep.py 重建 network.json。
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

HERE = Path(__file__).parent
DATASET_FILE = HERE / "tcmsp_dataset.json"

ALIASES = {
    "herb": ["herb", "中药", "中药名", "药材", "药材名", "drug", "herb_name"],
    "mol_id": ["mol_id", "molid", "成分id", "成分编号", "molecule_id", "id"],
    "name": ["name", "成分", "成分名", "化合物", "molecule_name", "molecule"],
    "en": ["en", "英文名", "english_name", "english", "mol_name"],
    "ob": ["ob", "口服生物利用度", "bioavailability", "oral_bioavailability"],
    "dl": ["dl", "类药性", "drug_likeness", "druglikeness"],
    "targets": ["targets", "靶点", "靶点基因", "gene", "genes", "target", "symbol"],
    "smiles": ["smiles", "canonical_smiles", "结构式"],
}


def match_column(columns: List[str], key: str) -> Optional[str]:
    lowered = {str(c).strip().lower(): c for c in columns}
    for alias in ALIASES[key]:
        if alias.lower() in lowered:
            return lowered[alias.lower()]
    return None


def split_targets(value) -> List[str]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).replace("，", ",").replace("、", ",")
    return [t.strip() for t in text.replace(";", ",").split(",") if t.strip()]


def load_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in (".xlsx", ".xls"):
        return pd.read_excel(path)
    for encoding in ("utf-8-sig", "utf-8", "gbk"):
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError:
            continue
    raise SystemExit(f"无法读取文件编码：{path}")


def build_dataset(df: pd.DataFrame) -> Dict[str, dict]:
    cols = {key: match_column(list(df.columns), key) for key in ALIASES}
    if not cols["herb"] or not cols["name"]:
        raise SystemExit(f"缺少必要列（中药名 / 成分名），识别到的列：{cols}")

    compounds: Dict[str, dict] = {}
    herb_compounds: Dict[str, List[str]] = {}

    for _, row in df.iterrows():
        herb = str(row[cols["herb"]]).strip()
        name = str(row[cols["name"]]).strip()
        if not herb or not name or herb.lower() == "nan":
            continue

        mol_id = str(row[cols["mol_id"]]).strip() if cols["mol_id"] else ""
        if not mol_id or mol_id.lower() == "nan":
            mol_id = f"IMP{abs(hash((herb, name))) % 10**8:08d}"

        compound = {
            "name": name,
            "en": str(row[cols["en"]]).strip() if cols["en"] else "",
            "ob": float(row[cols["ob"]]) if cols["ob"] and pd.notna(row[cols["ob"]]) else 0.0,
            "dl": float(row[cols["dl"]]) if cols["dl"] and pd.notna(row[cols["dl"]]) else 0.0,
        }
        if cols["smiles"] and pd.notna(row[cols["smiles"]]):
            compound["smiles"] = str(row[cols["smiles"]]).strip()
        if cols["targets"]:
            targets = split_targets(row[cols["targets"]])
            if targets:
                compound["targets"] = targets

        compounds[mol_id] = compound
        herb_compounds.setdefault(herb, [])
        if mol_id not in herb_compounds[herb]:
            herb_compounds[herb].append(mol_id)

    return {"COMPOUNDS": compounds, "HERB_COMPOUNDS": herb_compounds}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="TCMSP 导出的 CSV / Excel 文件")
    parser.add_argument("--replace", action="store_true",
                        help="只使用导入的数据，不使用内置示例数据")
    args = parser.parse_args()

    path = Path(args.file)
    if not path.exists():
        raise SystemExit(f"文件不存在：{path}")

    df = load_table(path)
    dataset = build_dataset(df)
    dataset["replace_builtin"] = args.replace
    DATASET_FILE.write_text(json.dumps(dataset, ensure_ascii=False, indent=2),
                            encoding="utf-8")

    print(f"已写入 {DATASET_FILE}")
    print(f"中药 {len(dataset['HERB_COMPOUNDS'])} 味，成分 {len(dataset['COMPOUNDS'])} 个")
    print("下一步：python data_prep.py  然后  streamlit run app.py")


if __name__ == "__main__":
    main()
