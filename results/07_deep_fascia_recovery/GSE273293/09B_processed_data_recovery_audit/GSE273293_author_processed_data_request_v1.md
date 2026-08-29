# GSE273293 / PRJNA1206333 processed-data request template

Status: prepared only; **no email has been sent**.

Public corresponding-author addresses reported with the article:

- Liangjun Li: `Liliangjun1212@sina.com`
- Danling Wang: `danlingwang@usc.edu.cn`

## Recommended English email

**Subject:** Request for sample-resolved processed scRNA-seq data for GSE273293 / PRJNA1206333

Dear Professor Li and Professor Wang,

I am conducting a non-commercial reanalysis of the human deep-fascia single-cell dataset reported in your article, “Identification of pro-fibrotic cellular subpopulations in fascia of gluteal muscle contracture using single-cell RNA sequencing” (J Transl Med, 2025; DOI: 10.1186/s12967-024-05889-y).

The public GSE273293 processed archive currently appears to contain one matrix triplet labelled `GSM8425828_GMC01`, whereas PRJNA1206333 contains 14 raw sequencing records corresponding to 10 GMC and 4 control samples. To reproduce sample-level analyses without treating cells as independent biological replicates, could you please share any one of the following, if available?

1. A sample-resolved Seurat, SingleCellExperiment or AnnData object;
2. Per-sample processed count matrices for GMC01–GMC10 and Con01–Con04;
3. A table mapping each cell barcode in the processed matrix to its sample ID and condition;
4. Cell-type/subcluster annotations and the preprocessing/QC scripts used in the article.

A repository or temporary institutional download link would be sufficient. The files will be used only for academic, non-commercial reproducibility analysis, with full citation of your article and data accessions.

Thank you for your time and for making the sequencing records publicly available.

Sincerely,

[Your name]
[Institution, if applicable]
[Contact email]

## 中文参考稿

**主题：申请 GSE273293 / PRJNA1206333 样本分辨的单细胞处理数据**

李老师、王老师您好：

我正在开展一项非商业性的学术复现分析，使用贵团队发表于 Journal of Translational Medicine 的人类深筋膜单细胞研究数据（DOI：10.1186/s12967-024-05889-y）。

目前 GSE273293 的公开处理数据似乎只包含一个标记为 `GSM8425828_GMC01` 的矩阵三件套，而 PRJNA1206333 中有对应 GMC 10 例、对照 4 例的 14 条原始测序记录。为了以样本/供体为统计单位复现分析，避免把单个细胞误作独立生物学重复，想请问贵团队是否可以提供以下任一资料：

1. 保留样本身份的 Seurat、SingleCellExperiment 或 AnnData 对象；
2. GMC01–GMC10 和 Con01–Con04 的逐样本处理后计数矩阵；
3. 处理矩阵中 `cell barcode → sample ID/condition` 的对应表；
4. 论文使用的细胞类型/亚群注释以及预处理、质控代码。

如有公开仓库或临时机构下载链接也完全可以。相关文件仅用于非商业性学术复现，后续将完整引用论文和数据登录号。

感谢您的时间，也感谢贵团队公开原始测序数据。

此致
敬礼

[姓名]
[单位（如适用）]
[联系邮箱]

