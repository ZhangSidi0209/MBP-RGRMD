#!/bin/bash
set -euo pipefail

# ======================================================
# 数据库更新 + CheckM2 质控 + dRep 去冗余 + GTDB 注释 + GEM索引构建(95/99)
# ======================================================

# ===== 基础路径设置 =====
BASE_DIR="/mnt/project/zsd/pipeline"
ALL_MAG_DIR="${BASE_DIR}/all_MAGs"
NEW_MAG_DIR="${BASE_DIR}/new_MAGs"
CHECKM_DIR="${BASE_DIR}/checkm"
CHECKM_NEW_DIR="${BASE_DIR}/checkm_new"
EUK_DIR="${BASE_DIR}/Eukaryote"
NEW_EUK_DIR="${BASE_DIR}/new_Eukaryote"   # <- 去掉空格，避免路径问题
DB_DIR="/mnt/data/database"
REFMAG_DIR="${DB_DIR}/refmag"
DATE_SUFFIX="$(date +%Y%m%d)"
THREADS=192

mkdir -p "${ALL_MAG_DIR}" "${NEW_MAG_DIR}" "${CHECKM_DIR}" "${EUK_DIR}" "${NEW_EUK_DIR}" "${REFMAG_DIR}"

log() { echo "[$(date '+%F %T')] $*"; }

# ===== 小工具函数 =====
# 把 FASTA 重命名为 >文件名|原头，并把每条序列写回（保留多条序列）
# 若某条序列 > 30,000,000 bp，则切成多段，避免 GEM 32Mb 限制
reheader_and_append() {
  local in_fa="$1"
  local out_fa="$2"
  local prefix
  prefix="$(basename "$in_fa")"
  awk -v prefix="$prefix" -v maxlen=30000000 '
    BEGIN{hdr=""; seq=""}
    function flush_chunk(id, chunk, part,    i, L, line) {
      if (chunk == "") return
      printf(">%s|%s_part%d\n", prefix, id, part) >> out
      L = length(chunk)
      for (i=1; i<=L; i+=60) {
        line = substr(chunk, i, 60)
        print line >> out
      }
    }
    function flush_seq(id, s,    L, part, start, end) {
      if (s == "") return
      L = length(s)
      if (L <= maxlen) {
        printf(">%s|%s\n", prefix, id) >> out
        for (i=1; i<=L; i+=60) print substr(s, i, 60) >> out
      } else {
        part = 1
        for (start=1; start<=L; start+=maxlen) {
          end = start + maxlen - 1
          if (end > L) end = L
          flush_chunk(id, substr(s, start, end-start+1), part)
          part++
        }
      }
    }
    BEGINFILE{ out="'"$out_fa"'" }
    /^>/ {
      if (hdr != "") { flush_seq(hdr, seq) }
      hdr = substr($0,2); seq=""
      next
    }
    { gsub(/\r/,""); seq = seq $0 }
    END { if (hdr != "") { flush_seq(hdr, seq) } }
  ' "$in_fa"
}

# ===== Step 1. 检查新的 MAG 文件 =====
log "Step 1. 检查新的 MAG 文件..."
shopt -s nullglob
NEW_FILES=( "${NEW_MAG_DIR}"/*.fa )
if [ ${#NEW_FILES[@]} -eq 0 ]; then
  log "未检测到新的 MAG 文件。"
else
  log "检测到 ${#NEW_FILES[@]} 个新 MAG 文件。"
  # Step 2. 移动并重命名新 MAG
  log "Step 2. 移动并重命名新 MAG..."
  for file in "${NEW_FILES[@]}"; do
    base="$(basename "$file" .fa)"
    new_name="${base}_${DATE_SUFFIX}.fa"
    # 如有同名则加时分秒
    if [ -f "${ALL_MAG_DIR}/${new_name}" ]; then
      new_name="${base}_${DATE_SUFFIX}_$(date +%H%M%S).fa"
    fi
    mv -f "$file" "${ALL_MAG_DIR}/${new_name}"
    log "Moved: $file -> ${ALL_MAG_DIR}/${new_name}"
  done
fi

# ===== Step 3. CheckM2 质控 =====
log "Step 3. 运行 CheckM2 质控..."
# CheckM2 要求 CHECKM2DB 指向数据库目录（包含多个文件），而不是单个 .dmnd
export CHECKM2DB="${DB_DIR}/CheckM2_database"
mkdir -p "${CHECKM_NEW_DIR}"

# 如果 ALL_MAG_DIR 中没有 .fa，则直接跳过，避免报错
if compgen -G "${ALL_MAG_DIR}/*.fa" > /dev/null; then
  conda run -n checkm2 checkm2 predict \
    --threads "${THREADS}" \
    --input "${ALL_MAG_DIR}/" \
    --output-directory "${CHECKM_NEW_DIR}" \
    -x fa --force
else
  log "ALL_MAG_DIR 中无 fasta，跳过 CheckM2。"
fi

# ===== Step 4. 提取高质量 MAG 名单 =====
log "Step 4. 筛选高质量 MAG..."
input_file_new="${CHECKM_NEW_DIR}/quality_report.tsv"
output_file_new="${CHECKM_NEW_DIR}/newlist_new"
if [ -f "${input_file_new}" ]; then
  awk 'BEGIN{FS="\t"} NR==1{next} {if($2>=80 && $3<=10) print $1}' "${input_file_new}" > "${output_file_new}"
  HQ_COUNT=$(wc -l < "${output_file_new}" || echo 0)
  log "高质量 MAG 数量：${HQ_COUNT}"
else
  log "未找到 ${input_file_new}，可能未运行或无新 MAG。创建空名单。"
  : > "${output_file_new}"
fi

# ===== Step 5. 复制高质量 MAG =====
log "Step 5. 复制高质量 MAG..."
mkdir -p "${CHECKM_DIR}/combined_filteredfa"
if [ -s "${output_file_new}" ]; then
  while IFS= read -r line; do
    src="${ALL_MAG_DIR}/${line}.fa"
    dst="${CHECKM_DIR}/combined_filteredfa/${line}.fa"
    if [ -f "$src" ]; then
      cp -f "$src" "$dst"
    else
      log "Warning: ${line}.fa not found in all_MAGs" >&2
    fi
  done < "${output_file_new}"
else
  log "高质量名单为空，跳过复制。"
fi

# ===== Step 6. dRep 去冗余（99/95）=====
log "Step 6. 运行 dRep 去冗余..."
for ANI in 99 95; do
  DREP_OUT="${BASE_DIR}/drep${ANI}"
  mkdir -p "${DREP_OUT}"
  GENOME_LIST="${CHECKM_DIR}/combined_filteredfa/genome_list_${ANI}.txt"
  find "${CHECKM_DIR}/combined_filteredfa" -type f -name "*.fa" > "${GENOME_LIST}"

  if [ ! -s "${GENOME_LIST}" ]; then
    log "GENOME_LIST 为空，跳过 dRep ${ANI}。"
    continue
  fi

  # 使用清单文件避免 “Argument list too long”
  dRep dereplicate "${DREP_OUT}" \
    --genomes "${GENOME_LIST}" \
    -p "${THREADS}" \
    --ignoreGenomeQuality \
    -sa 0.${ANI} || true  # dRep 某些绘图报错不影响 derep 结果，容忍返回码
done

# ===== Step 7. GTDB-Tk 分类 =====
log "Step 7. 运行 GTDB-Tk 分类..."
export GTDBTK_DATA_PATH="${DB_DIR}/gtdb/release220/"
for ANI in 99 95; do
  DREP_DIR="${BASE_DIR}/drep${ANI}/dereplicated_genomes"
  if compgen -G "${DREP_DIR}/*.fa" > /dev/null; then
    gtdbtk classify_wf \
      --skip_ani_screen \
      --genome_dir "${DREP_DIR}" \
      --out_dir "${BASE_DIR}/gtdb${ANI}" \
      -x fa --cpus "${THREADS}"
  else
    log "未找到 ${DREP_DIR}/*.fa，跳过 GTDB ${ANI}。"
  fi
done

# ===== Step 8. 整理数据库文件 =====
log "Step 8. 更新数据库文件..."
mkdir -p "${CHECKM_DIR}/diamond_output" "${CHECKM_DIR}/protein_files"

if [ -d "${CHECKM_NEW_DIR}/diamond_output" ]; then
  shopt -s nullglob
  for file in "${CHECKM_NEW_DIR}/diamond_output"/*; do
    base="$(basename "$file")"
    if [[ -f "${CHECKM_DIR}/diamond_output/$base" ]]; then
      mv -f "$file" "${CHECKM_DIR}/diamond_output/${base%.*}_new.${base##*.}"
    else
      mv -f "$file" "${CHECKM_DIR}/diamond_output/"
    fi
  done
fi

if [ -d "${CHECKM_NEW_DIR}/protein_files" ]; then
  shopt -s nullglob
  mv -f "${CHECKM_NEW_DIR}/protein_files"/* "${CHECKM_DIR}/protein_files/" 2>/dev/null || true
fi

[ -f "${CHECKM_NEW_DIR}/quality_report.tsv" ] && cp -f "${CHECKM_NEW_DIR}/quality_report.tsv" "${CHECKM_DIR}/quality_report.tsv"
[ -f "${CHECKM_NEW_DIR}/newlist_new" ] && cp -f "${CHECKM_NEW_DIR}/newlist_new" "${CHECKM_DIR}/newlist"

rm -rf "${CHECKM_NEW_DIR}" || true
log "数据库更新完成，新 MAG 已整合。"

# ===== Step 9. 真核 MAG 更新 =====
log "Step 9. 更新真核 MAG..."
if compgen -G "${NEW_EUK_DIR}/*.fa" > /dev/null; then
  for file in "${NEW_EUK_DIR}"/*.fa; do
    base="$(basename "$file" .fa)"
    new_name="${base}_${DATE_SUFFIX}.fa"
    if [ -f "${EUK_DIR}/${new_name}" ]; then
      new_name="${base}_${DATE_SUFFIX}_$(date +%H%M%S).fa"
    fi
    mv -f "$file" "${EUK_DIR}/${new_name}"
    log "Moved new Euk: $file -> ${EUK_DIR}/${new_name}"
  done
else
  log "未检测到新的真核 MAG。"
fi

# ===== Step 10. 构建 GEM 索引（95% 和 99%）=====
log "Step 10. 构建 GEM 索引..."
for ANI in 95 99; do
  REF_DIR="${BASE_DIR}/ref${ANI}"
  TMP_DIR="${BASE_DIR}/tmp_ref${ANI}_${DATE_SUFFIX}"
  MERGED_FA="${REF_DIR}/ref${ANI}.nr.fa"
  mkdir -p "${REF_DIR}" "${TMP_DIR}"

  # 复制去冗余 MAG 到 tmp
  if compgen -G "${BASE_DIR}/drep${ANI}/dereplicated_genomes/*.fa" > /dev/null; then
    cp -f "${BASE_DIR}/drep${ANI}/dereplicated_genomes/"*.fa "${TMP_DIR}/"
  else
    log "警告：未找到 drep${ANI} 去冗余 MAG，继续构建（仅真核/历史数据）。"
  fi

  # 复制真核 MAG 到 tmp（若存在）
  if compgen -G "${EUK_DIR}/*.fa" > /dev/null; then
    cp -f "${EUK_DIR}/"*.fa "${TMP_DIR}/"
  fi

  # 合并并重命名（避免 >32Mb 单序列；必要时分块）
  : > "${MERGED_FA}"
  shopt -s nullglob
  for f in "${TMP_DIR}"/*.fa; do
    log "合并引用：$(basename "$f")"
    reheader_and_append "$f" "${MERGED_FA}"
  done

  # 构建 GEM 索引
  mkdir -p "${REFMAG_DIR}"
  gem-indexer -i "${MERGED_FA}" -o "${REFMAG_DIR}/ref${ANI}" -t "${THREADS}"
  log "[DONE] GEM ${ANI}% 索引生成: ${REFMAG_DIR}/ref${ANI}.gem"

  # 清理临时目录（保留 MERGED_FA 便于溯源）
  rm -rf "${TMP_DIR}"
done

log "全部流程完成"
