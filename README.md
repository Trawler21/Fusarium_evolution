# Fusarium_evolution

本仓库包含从 OrthoFinder 单拷贝 orthogroup 提取序列并构建系统发育树的脚本。

## 原有脚本
- `make_sc_og_species_gene_table.py`：从 `Orthogroups.tsv` 和 `Orthogroups_SingleCopyOrthologues.txt` 生成单拷贝 OG 的宽表/长表。
- `extract_sc_ogs_by_table.py`：根据长表从各物种蛋白组中提取 OG 序列（支持覆盖度过滤）。
- `check_og_coverage_and_dup.py`：检查每个 OG 是否存在缺失/重复/额外物种。
- `check_trimmed_aln_lengths.py`：检查每个裁剪后比对文件内部序列长度是否一致。
- `concat_trimmed_alignments.py`：串联所有裁剪后的比对并输出分区信息（已按你上传到 GitHub 的脚本作为当前流程使用版本）。

## Nextflow 流程（新增）
新增 `main.nf + nextflow.config`，将完整流程封装为可复现的 pipeline：

1. 生成单拷贝表（wide + long）
2. 提取 OG 蛋白 FASTA
3. 检查 OG 覆盖度与重复
4. MAFFT 比对
5. trimAl 裁剪
6. 裁剪后长度一致性检查
7. 串联比对
8. （可选）IQ-TREE 建树

### 依赖
- Python 3
- Python 包：`pandas`、`biopython`
- `mafft`
- `trimal`
- `nextflow`
- （可选）`iqtree2`

### 运行示例
```bash
nextflow run main.nf \
  --orthogroups_tsv /path/to/Orthogroups.tsv \
  --single_copy_txt /path/to/Orthogroups_SingleCopyOrthologues.txt \
  --proteomes_dir /path/to/proteomes \
  --species_ref /path/to/species_193.txt \
  --expect_species_n 193 \
  --total_species 193 \
  --min_coverage 0.5 \
  --outdir results \
  --do_iqtree false \
  --concat_script /path/to/concat_trimmed_alignments.py
```

开启 IQ-TREE：
```bash
nextflow run main.nf \
  --orthogroups_tsv /path/to/Orthogroups.tsv \
  --single_copy_txt /path/to/Orthogroups_SingleCopyOrthologues.txt \
  --proteomes_dir /path/to/proteomes \
  --species_ref /path/to/species_193.txt \
  --do_iqtree true \
  --iqtree_extra "-T AUTO"
```


- `--concat_script`：指定串联脚本路径（默认使用仓库内 `concat_trimmed_alignments.py`）。

### 输出目录
- `results/aln/`：MAFFT 输出
- `results/trimmed/`：trimAl 输出
- `results/concat/concatenated_alignment.fasta`：串联矩阵
- `results/concat/partition.tsv`：位点分区
- `results/iqtree/`：IQ-TREE 输出（若开启）
- 以及 Nextflow trace/report/timeline/dag

## 你还需要确认的信息
为了在你的数据上直接跑通，建议先确认以下参数：
- 物种蛋白 FASTA 文件命名是否是 `物种名.xxx`（脚本会以第一个 `.` 前字符串作为物种名）。
- `species_193.txt` 中的物种名是否与 FASTA 文件前缀完全一致。
- 是否要在串联前过滤掉长度过短的 OG（当前流程不做此过滤）。
- IQ-TREE 的模型搜索与 bootstrap 策略（当前默认 `-m MFP -bb 1000 -alrt 1000`）。
