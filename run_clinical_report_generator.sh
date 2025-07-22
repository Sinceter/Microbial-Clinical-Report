#!/usr/bin/env bash
set -e

#
# run_clinical_report_generator.sh
#
# 假设目录结构：
#   Clinical_Report_Generator/
#   ├── full_report_template.html
#   ├── generate_full_report.py
#   ├── generate_sliding_tables.py
#   ├── mNGS_Kraken-based.R
#   ├── sliding_tables_template.html
#   ├── run_clinical_report_generator.sh   ← （本脚本）
#   └── Functions/                          ← R 脚本内部 source 这些辅助脚本
#       ├── Extract_Nanopore_Info.R
#       ├── Handle_Kraken_Style_Reports.R
#       ├── Match_Judge.R
#       └── Naming_Correction.R
#
# 用法示例：
#   cd Clinical_Report_Generator
#   ./run_clinical_report_generator.sh \
#     --filepath ./test/mixed_nanopore_run \
#     --culture_table ./clinical_culture.tsv \
#     --col_header_NR RunID \
#     --col_header_Bc Barcode \
#     --col_header_Cl ClinicalResult \
#     --read_count_cutoff 1 \
#     --read_perct_cutoff 10,20 \
#     --show_top 50 \
#     --time 12 \
#     --output ./my_output_dir \
#     --outname MyCustomReportName
#
# 关键说明：
#   1. 如果用户传了 "--output <dir>"，则直接使用该目录，否则，默认输出到
#        dirname(filepath)/output_<YYYYMMDD-HHMM>
#   2. 最终会在该目录下寻找：
#        Table_Clean_mNGS_<basename(filepath)>_<YYYYMMDD>.tsv
#   3. 其余逻辑同之前版本；支持 --outname/-n 覆盖最终报告基名。
#

# -------- 帮助信息 --------
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
  sed -n '1,50p' "$0"
  exit 0
fi

# 保存所有参数
ARGS=("$@")

# -------- 解析 --outname / -n 参数 --------
OUTNAME=""
REMAINING_ARGS=()

i=0
while [[ $i -lt $# ]]; do
  case "${ARGS[i]}" in
    --outname|-n)
      OUTNAME="${ARGS[i+1]}"
      i=$((i+2))
      ;;
    *)
      REMAINING_ARGS+=("${ARGS[i]}")
      i=$((i+1))
      ;;
  esac
done

# -------- 步骤 1：运行 R 脚本 mNGS_Kraken-based.R --------
echo ">> Running R script: mNGS_Kraken-based.R ${REMAINING_ARGS[*]}"
Rscript mNGS_Kraken-based.R "${REMAINING_ARGS[@]}"

# -------- 步骤 2：提取 filepath 和 output_dir --------

# 从 REMAINING_ARGS 中提取 --filepath 或 -f
INPUT_PATH_RAW=""
for ((j=0; j<${#REMAINING_ARGS[@]}; j++)); do
  if [[ "${REMAINING_ARGS[j]}" == "--filepath" || "${REMAINING_ARGS[j]}" == "-f" ]]; then
    INPUT_PATH_RAW="${REMAINING_ARGS[j+1]}"
    break
  fi
done

if [[ -z "$INPUT_PATH_RAW" ]]; then
  echo "Error: --filepath (or -f) must be provided."
  exit 1
fi

# 去掉末尾可能的斜杠
INPUT_PATH="${INPUT_PATH_RAW%/}"

# 从 REMAINING_ARGS 中提取 --output 或 -o（可选）
OUTPUT_DIR=""
for ((j=0; j<${#REMAINING_ARGS[@]}; j++)); do
  if [[ "${REMAINING_ARGS[j]}" == "--output" || "${REMAINING_ARGS[j]}" == "-o" ]]; then
    OUTPUT_DIR="${REMAINING_ARGS[j+1]}"
    break
  fi
done

# 如果没有显式提供 "--output"，则默认使用 dirname(filepath)/output_<YYYYMMDD-HHMM>
if [[ -z "$OUTPUT_DIR" ]]; then
  PARENT_DIR=$(dirname "$INPUT_PATH")
  TIMESTAMP=$(date +%Y%m%d-%H%M)
  OUTPUT_DIR="${PARENT_DIR}/output_${TIMESTAMP}"
fi

# 确保 OUTPUT_DIR 存在
if [[ ! -d "$OUTPUT_DIR" ]]; then
  echo "Error: Cannot find R script output directory: $OUTPUT_DIR"
  exit 1
fi
echo ">> R script output directory: $OUTPUT_DIR"

# -------- 步骤 3：定位 TSV 文件 --------

BASE=$(basename "$INPUT_PATH")
DATE=$(date +%Y%m%d)
TSV_NAME="Table_Clean_mNGS_${BASE}_${DATE}.tsv"
TSV_PATH="${OUTPUT_DIR}/${TSV_NAME}"

if [[ ! -f "$TSV_PATH" ]]; then
  echo "Error: Cannot find TSV at expected path: $TSV_PATH"
  exit 1
fi
echo ">> Found TSV: $TSV_PATH"

# -------- 步骤 4：生成临时 sliding tables HTML 文件 --------

TMP_BASE=$(mktemp -u /tmp/sliding_tables_XXXXXX)
TMP_TABLES="${TMP_BASE}.html"
echo ">> Temporary sliding tables file: $TMP_TABLES"

# -------- 步骤 5：调用 generate_sliding_tables.py --------

echo ">> Generating sliding tables..."
python3 generate_sliding_tables.py \
  --tsv "$TSV_PATH" \
  --template sliding_tables_template.html \
  --output "$TMP_TABLES"

if [[ ! -f "$TMP_TABLES" ]]; then
  echo "Error: generate_sliding_tables.py did not create $TMP_TABLES"
  exit 1
fi
echo ">> Sliding tables saved to: $TMP_TABLES"

# -------- 步骤 6：确定最终报告文件名 --------

if [[ -n "$OUTNAME" ]]; then
  FINAL_BASE="$OUTNAME"
else
  TS2=$(date +%Y%m%d%H%M)
  FINAL_BASE="final_clinical_report_${TS2}"
fi
FINAL_HTML="${FINAL_BASE}.html"

# -------- 步骤 7：调用 generate_full_report.py --------

echo ">> Generating final HTML report: $FINAL_HTML"
python3 generate_full_report.py \
  --tables "$TMP_TABLES" \
  --template full_report_template.html \
  --output "$FINAL_HTML"

if [[ ! -f "$FINAL_HTML" ]]; then
  echo "Error: generate_full_report.py did not create $FINAL_HTML"
  rm -f "$TMP_TABLES"
  exit 1
fi
echo ">> Final HTML created: $FINAL_HTML"

# -------- 步骤 8：删除临时滑动表格文件 --------

rm -f "$TMP_TABLES"
echo ">> Removed temporary file: $TMP_TABLES"

# -------- 步骤 9：使用 wkhtmltopdf & wkhtmltoimage 生成 PDF 和 PNG --------

PDF_OUT="${FINAL_BASE}.pdf"
PNG_OUT="${FINAL_BASE}.png"

echo ">> Converting HTML to PDF: $PDF_OUT"
wkhtmltopdf "$FINAL_HTML" "$PDF_OUT"

echo ">> Converting HTML to PNG: $PNG_OUT"
wkhtmltoimage "$FINAL_HTML" "$PNG_OUT"

echo "✔ Pipeline complete!"
echo "   • HTML: $FINAL_HTML"
echo "   • PDF : $PDF_OUT"
echo "   • PNG : $PNG_OUT"
