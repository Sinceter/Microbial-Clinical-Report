#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
generate_sliding_tables.py

功能：
  - 读取一个列名固定的 TSV 文件（含 “Nanopore run, Time, Barcode, Database, CS_score, Software, ...” 等列）。
  - 按 Software + Database + CS_score 的组合分组。
  - 每个组合生成一个“滑动表格”片段，标题为：
        <Software> Report (Database: <Database>, CS_score: <CS_score>)
  - 输出到一个 HTML 文件（只包含所有滑动表格片段）。

使用示例：
  python3 generate_sliding_tables.py \
      --tsv Table_Clean_mNGS.tsv \
      --template sliding_tables_template.html \
      --output tables_output.html
"""

import argparse
import pandas as pd
from jinja2 import Template


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate multiple sliding tables (grouped by Software, Database, CS_score) from a TSV"
    )
    parser.add_argument(
        "--tsv", "-i", required=True,
        help="Path to the input TSV file (column names must match exactly)"
    )
    parser.add_argument(
        "--template", "-t", required=True,
        help="Path to the sliding table Jinja2 template (sliding_tables_template.html)"
    )
    parser.add_argument(
        "--output", "-o", required=True,
        help="Path to the output HTML file (e.g. tables_output.html)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # 1. 读取 TSV，所有列都当字符串，缺失填 ""
    df = pd.read_csv(args.tsv, sep="\t", dtype=str).fillna("")

    # 检查关键列是否存在
    required_cols = ["Software", "Database", "CS_score"]
    for col in required_cols:
        if col not in df.columns:
            print(f"ERROR: Column '{col}' not found in TSV. Please check column names.")
            return

    # 记录所有列名（保留顺序）
    columns = list(df.columns)

    # 2. 按 [Software, Database, CS_score] 分组
    grouped = df.groupby(["Software", "Database", "CS_score"], sort=False)

    sections = []
    for (software_name, database_val, cs_score_val), subdf in grouped:
        # subdf 是该组合下的所有行
        # 转成记录行字典列表
        rows = subdf.to_dict(orient="records")

        sections.append({
            "software": software_name,
            "database": database_val if database_val.strip() != "" else None,
            "cs_score": cs_score_val if cs_score_val.strip() != "" else None,
            "rows": rows
        })

    # 3. 渲染 Jinja2 模板
    with open(args.template, "r", encoding="utf-8") as fp:
        tpl = Template(fp.read())

    rendered = tpl.render(
        columns=columns,
        sections=sections
    )

    # 4. 输出到 HTML
    with open(args.output, "w", encoding="utf-8") as fp:
        fp.write(rendered)

    print(f"✔ Sliding tables generated: {args.output}")


if __name__ == "__main__":
    main()
