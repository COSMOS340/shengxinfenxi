# ACTIVE - 停止并清理 PRJNA932556，切换到 GSE205506

状态：`REQUESTED`

日期：2026-07-16

## 目标

立即停止 `PRJNA932556` 的 Cell Ranger 主线，删除该队列占用硬盘的大体积数据，然后使用 `GSE205506` 作者处理矩阵建立 CRC 正式 R/Seurat 对象。

## 第一阶段：停止和删除 PRJNA932556

1. 查找并停止所有正在读取或写入 `PRJNA932556`、`SRR23490337`、`R1-R3` 或 `S1-S3` 数据的 Cell Ranger、QEMU、Linux 子系统或辅助进程。
2. 不得根据名称相似性直接删除。先从现有请求、运行日志和实际进程命令行提取精确绝对路径，生成 `prjna932556_deletion_manifest_before.tsv`。
3. 清单字段必须包含：`absolute_path`、`path_type`、`size_bytes`、`evidence_linking_path_to_prjna932556`、`deletion_status`、`error_message`。
4. 删除确认属于 `PRJNA932556` 的 BAM、BAI、FASTQ、Cell Ranger 工作目录、临时文件和表达输出。
5. 保留脚本、manifest、校验和、运行日志和轻量状态文件。
6. 不得删除 Cell Ranger 软件、参考基因组、`GSE205506`、`GSE189926`、`GSE235863`、其他项目或无法明确归属的路径。
7. 删除后重新扫描删除前清单中的每个精确路径，生成 `prjna932556_deletion_audit_after.tsv`，并报告 `total_bytes_before`、`total_bytes_deleted`、`total_bytes_remaining`。
8. 对仍在占用、删除失败或归属不明确的路径逐项说明，不得使用通配符扩大删除范围。

## 第二阶段：获取 GSE205506

1. 从 GEO Series `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205506` 获取页面实际列出的 `GSE205506_RAW.tar`。
2. 记录最终下载地址、实际文件名、字节数和 SHA-256，不使用手工推导的下载地址替代实际地址。
3. 大文件放在 F 盘本项目的 `GSE205506` 专用目录。输出 `gse205506_input_inventory.tsv`，记录所有实际绝对路径。
4. 解包前检查可用空间。解包后导出 tar 内全部文件名、大小、矩阵格式和样本标识，不猜测文件命名规则。
5. 大矩阵和 RDS 保留在 PC，不通过 GitHub 上传。

## 第三阶段：精确疗效映射

1. 从 DOI `10.1016/j.ccell.2023.04.011` 的论文页面获取实际补充 `Table S1`。
2. 记录实际补充文件名、下载地址、字节数和 SHA-256。
3. 从表中按实际字段提取19位患者的 `pCR/non-pCR` 信息，不根据患者编号、样本顺序、正文图形或15比4的总数自行分配。
4. 与现有 GEO metadata 的 `subject` 字段精确连接，导出 `gse205506_response_mapping.tsv` 和 `gse205506_response_join_audit.tsv`。
5. 报告未匹配、重复、多对多连接和字段空值。存在任何未解释问题时，停止疗效分析并回传审计表。

## 第四阶段：R/Seurat 正式对象

1. 全程使用 R。禁止调用 Python、Scanpy或继承其他对象的 UMAP、cluster、双细胞和注释。
2. 从实际作者处理矩阵建立每样本 Seurat 对象，保留原始 GEO 字段和补充表字段，不覆盖源字段。
3. 在逐患者疗效映射通过前，只运行矩阵结构审计和不依赖疗效的基础 QC。
4. 映射通过后完成 QC、环境 RNA 审计、双细胞审计、标准化、PCA、批次结构评估、UMAP、多个聚类分辨率和 markers。
5. 导出每个 cluster 的 Top 10 和 Top 50 markers、注释依据、患者级细胞数和样本保留率。
6. QC 阈值依据每个样本实际分布确定并输出敏感性审计，不机械复用其他队列阈值。

## 必须回传到 GitHub 的轻量文件

- `STATUS.md`
- `prjna932556_deletion_manifest_before.tsv`
- `prjna932556_deletion_audit_after.tsv`
- `gse205506_input_inventory.tsv`
- `gse205506_tar_inventory.tsv`
- `gse205506_response_mapping.tsv`
- `gse205506_response_join_audit.tsv`
- `gse205506_matrix_structure_audit.tsv`
- `gse205506_qc_summary.tsv`
- `gse205506_cluster_top10_markers.tsv`
- `gse205506_cluster_top50_markers.tsv`
- `gse205506_annotation_evidence.tsv`
- R 脚本、参数文件、`sessionInfo()` 和 `output_manifest.tsv`
- 低体积 PNG/PDF QC 图、UMAP 和 marker dotplot

上传目录：

`bioinformatics-reproduction/projects/pan-gi-ici-ccc-20260709/desktop_exchange/uploads/20260716_stop_prjna932556_and_build_gse205506/`

## 验收条件

- `PRJNA932556` 大体积数据已按删除前清单逐项核验删除，释放空间有字节级记录。
- Cell Ranger 软件、参考基因组和其他数据集未被删除。
- `GSE205506` 来源、实际矩阵结构和逐患者疗效映射可追溯。
- 未核对的疗效字段没有进入任何图或统计。
- 正式对象和图均由本轮 R 流程生成。
