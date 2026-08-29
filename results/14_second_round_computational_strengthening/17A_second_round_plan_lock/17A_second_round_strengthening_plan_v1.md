# 第二轮干实验补强方案（正式锁定版）

## 0. 文档状态

- 项目：人类筋膜–肌肉机械敏感性研究
- 版本：v1.0，正式执行版
- 锁定日期：2026-08-26
- 方案依据：`14A_第二轮补强方案.md`、Step 12B、Step 12C、Step 15F 及当前论文补强版
- 当前总体证据等级：**CAUTION**
- 当前解释类别：`ACTOMYOSIN_PERTURBATION_PLAUSIBILITY_WITH_PROLIFERATION_COMPETITION`
- 本轮定位：以最少的新增分析收窄三个已明确的不确定性，不通过新增候选、放宽阈值或堆叠数据集来制造证据升级
- 执行原则：先审计来源和设计，再冻结分析契约，最后读取表达结果；所有负结果和不可估计结果均保留

## 1. 锁定后的研究目标

本轮不再试图证明“普遍的单基因机械敏感性机制”，而是回答三个具体问题：

1. 机械相关程序在单细胞层面是否存在稳定的**共表达层级证据**，还是主要由细胞群组成差异造成？
2. 在无机械加载条件的外部转录组中，TGFβ 暴露效应与 TEAD 抑制效应能否支持机械相关程序的**调控轴交叉验证**？
3. 在独立的跨组织数据中，机械相关模块是否呈现刚度相关的**剂量–反应形式**，以及是否存在可重复的疾病遗传学关联？

本轮结果的价值是缩小解释边界、排序功能实验优先级和改善论文透明度；不替代 Step 16A 的供体重复机械功能实验。

## 2. 必须保留的核心结论与边界

### 2.1 可以保留

- GSE300230 在两个培养人真皮成纤维细胞背景中支持 actomyosin 机械张力相关转录响应的 perturbation plausibility。
- 人类成纤维细胞状态与多组分 ECM/integrin–cytoskeletal 程序存在情境依赖的关联。
- Integrin/focal-adhesion 与 Rho/actomyosin 是 Step 16A 中优先进行功能扰动的候选分支。

### 2.2 本轮不得直接宣称

- `cell-intrinsic mechanosensitivity` 已被证明；
- YAP/TAZ–TEAD 独立驱动了机械程序；
- GSE123100、GSE276045/46 或 GSE338388 是筋膜的直接复制；
- 纤维肌痛 GWAS 富集等同于疼痛因果机制；
- PIEZO2 是普遍或必需的机械感受器；
- S1–S3 全部阳性即可取消 CAUTION；
- 公共数据可以替代至少 3 个独立供体的机械功能扰动实验。

### 2.3 本轮证据升级上限

即使 S1–S3 均得到支持，论文最多使用：

> **evidence-bounded mechanistic hypothesis with computational triangulation**

整体证据等级仍保持 **CAUTION**，直到 Step 16A 满足供体重复、功能表型、通路扰动和混杂控制等预设条件。

## 3. 已纠正的数据集定位

这是本轮执行前的强制修正：

| 数据集 | 正式定位 | 不可作何种解释 |
|---|---|---|
| PRJNA607098 | S1 单细胞层共表达与样本内组成效应审计 | 不直接证明细胞内在因果机制 |
| GSE130973 | S1 的独立皮肤成纤维细胞补充复核 | 不与 PRJNA607098 跨研究合并成一个统计样本 |
| GSE338388 | S2：TGFβ 暴露 × TEAD 抑制的调控轴交叉验证 | 无机械加载/刚度条件，不能证明机械因果或独立机械驱动 |
| **GSE123100** | **S3 正式刚度剂量–反应主力：HTM 细胞刚度梯度** | 不能作为筋膜直接验证；只检验剂量–反应形式 |
| GSE276045（及其伴随 accession，如审计确认） | S3 的 WI-38 肺成纤维细胞刚度 × 增殖复核与混杂检查 | 不应简单写成“另一套刚度梯度主力”；不能作筋膜直接验证 |
| GCST90838603–11 | S3 可选遗传学正交关联 | 不得写作因果通路或细胞类型因果证明 |

## 4. Frozen 分析契约

### 4.1 全局规则

- 沿用原 15-gene frozen candidate panel 和既有 frozen modules；不新增候选，不替换缺失基因，不事后调整阈值。
- 各研究分别分析；不把不同研究、组织、平台的表达矩阵直接拼接后作单一模型。
- 细胞是嵌套观测，donor/sample 是生物学统计单位；不得把细胞数量当作独立重复数。
- 主结果报告效应量、方向、不确定性、覆盖率和缺失情况，不只报告 p 值。
- 模块基因不足时标记为 `NOT_ESTIMABLE` 或 `PARTIAL_COVERAGE`，不进行插补或候选替换。
- 预设分析失败不转化为“换数据继续寻找阳性”；失败本身进入证据边界和停止记录。

### 4.2 S1：单细胞共表达层级审计

**主要输入：** PRJNA607098 目标细胞/样本的现有远程 Zarr 流式通道；GSE130973 已有 Seurat 对象及其受试者字段。

**主要问题：** 在每个样本内，integrin/focal-adhesion 分支与 actomyosin/机械程序的共检测是否超出由测序深度、基因检出率和细胞组成所产生的随机背景？

**锁定做法：**

1. 先锁定目标细胞状态、样本字段、基因覆盖率和每个样本的有效细胞数。
2. 双阳性比例仅作为描述性细胞层指标；主要推断在样本/供体层汇总。
3. 采用样本内置换或条件化随机背景，至少控制 library size、gene detection rate 和目标基因的边际检出率。
4. 以样本为单位汇总共表达效应；必要时用带样本随机效应的二项/β-二项模型，避免伪重复。
5. 在样本内或 donor-level 层面进行竞争程序敏感性分析，至少保留 cell-cycle、ECM/TGF、炎症/缺氧等既有竞争解释。
6. PRJNA607098 与 GSE130973 分开输出；只做方向和证据层级比较，不做跨研究合并检验。

**S1 允许的结论：**

- `cell-level co-expression evidence`；
- `composition-consistent`、`mixed` 或 `within-sample co-expression supported` 等层级描述；
- 若覆盖或样本数不足，则为 `NOT_ESTIMABLE`。

**S1 禁止的结论：** `cell-intrinsic mechanism proven`、细胞内因果、筋膜特异性。

### 4.3 S2：TGFβ 暴露 × TEAD 抑制调控轴交叉验证

**前置条件：** GSE338388 必须先通过官方来源、样本设计、处理标签、重复结构和 frozen panel/module 覆盖率审计；审计未通过则不进入表达分析。

**锁定设计：** TGFβ 暴露（±）× TEAD 抑制（±）的双因素框架。主效应和交互项必须分别报告。

**锁定分析：**

- YAP/TAZ–TEAD 相关靶模块与 TGFβ–SMAD 相关靶模块分别定义并冻结；
- 先检验 TEAD 抑制对 YAP/TAZ 模块的特异性影响，再检验 TGFβ 暴露对 SMAD 模块的影响；
- frozen actomyosin、integrin/focal-adhesion、ECM、cell-cycle 和炎症模块作为结果模块并列报告；
- 将“轴特异性响应”与“机械程序归属”分开：即使某模块主要受 TEAD 抑制影响，也只能支持调控轴交叉验证；
- 报告跨组织来源、样本量、覆盖率和处理结构限制。

**S2 允许的结论：** `TGFβ/SMAD and TEAD-related regulatory-axis cross-validation`，以及在该数据集内的模块方向和交互模式。

**S2 禁止的结论：** 无机械条件下的机械因果、YAP/TAZ 独立驱动已被证明、筋膜直接复制。

### 4.4 S3：刚度剂量–反应与遗传学正交

#### S3.1 GSE123100：正式剂量–反应主力

- 先审计刚度单位、实际梯度、重复结构、批次和细胞状态；
- 以刚度的预设连续变换作为主自变量，优先检验 module score 与刚度的方向趋势；
- 非线性形态、单调性和候选基因结果作为次级描述；
- 报告模块覆盖、效应大小、置信区间/不确定性和样本级统计；
- 结果定位为跨组织的剂量–反应形式检验，不写成筋膜直接验证或因果机制。

#### S3.2 GSE276045（及经审计确认的伴随 accession）：刚度 × 增殖复核

- 先确认其确为 WI-38 人肺成纤维细胞刚度 × 增殖设计；
- 将增殖状态作为设计因素或预设协变量，检验刚度方向是否在增殖层级间一致；
- 若结构不支持独立估计，则标记 `NOT_ESTIMABLE`，不强行校正；
- 该资源用于独立跨组织复核和竞争解释，不取代 GSE123100 的正式剂量–反应主分析。

#### S3.3 GWAS：可选正交层

- 仅在 summary statistics、性状定义、版本和 LD 参考可追溯时执行；
- 预先固定 gene set、基因映射、LD 参考、基因集大小处理、基因间相关控制和多重检验方案；
- 优先使用 MAGMA competitive gene-set test；stratified LDSC 只作为条件允许时的补充；
- 逐个性状报告结果，不把多个 GWAS 条目挑选成单一阳性结论；
- 正结果只能支持疾病相关性锚点，不能支持纤维肌痛/疼痛因果通路。

## 5. 执行顺序与门控

### Gate 17A：方案与主张锁定（本步骤）

状态：**PASS，已锁定。** 已完成数据集主次纠正、主张边界、统计单位和停止规则固定。

### Gate 17B：S1 来源与可分析性审计

优先级最高，原因是零新增大文件获取，直接回应 Step 15D 的 `NOT_ESTIMABLE`。只有在目标细胞、样本字段、基因覆盖和有效样本数满足要求时，才进入 S1 分析。

### Gate 17C：S1 共表达层级分析

先完成描述性共检测，再完成样本内置换、检测率校正和 donor/sample-level 汇总。输出不得使用 `cell-intrinsic` 作为自动判定标签。

对应 R 脚本：`scripts/17C_analyze_S1_coexpression_level.R`；远程定向读取助手：`scripts/17C_extract_PRJNA607098_S1_coexpression.py`。

### Gate 17D：S2 GSE338388 审计及双因素分析

先 provenance/design audit，后 frozen panel/module 分析。若实际设计不是 TGFβ ± × TEAD 抑制 ±，则降级为描述性或停止。

对应审计脚本：`scripts/17D_audit_GSE338388_TGFb_TEAD_design.R`；人工复核后的标签校正脚本：`scripts/17D2_correct_GSE338388_design_labels.R`；矩阵准入审计脚本：`scripts/17D3_audit_GSE338388_processed_matrix.R`。

### Gate 17E：S3 GSE123100/GSE276045 审计及剂量–反应

GSE123100 为主，GSE276045/伴随 accession 为复核。GWAS 为可选分支，不作为核心计算补强能否完成的硬门槛。

### Gate 17F：证据重综合

统一更新证据域摘要、主张–证据映射、论文 Results/Discussion 和限制性声明。无论本轮阳性数量如何，整体证据等级默认为 CAUTION，除非未来供体重复功能实验另行满足 Step 12C/16A 预设条件。

## 6. 风险、停止与退化规则

| 风险/情形 | 处理规则 |
|---|---|
| PRJNA607098 无法稳定取得样本字段或目标细胞 | S1 标为 `NOT_ESTIMABLE`；不以细胞级伪重复替代 |
| 单细胞基因检出率不足 | 只报告覆盖率和描述性结果；不把缺失当阴性 |
| GSE130973 donor mapping 不完整 | 不做 donor-level 推断；保留为受限补充来源 |
| GSE338388 不是双因素设计或处理标签含糊 | 停止轴间推断，最多保留来源审计结果 |
| GSE123100 梯度或重复结构不支持剂量模型 | 标为 `NOT_ESTIMABLE`，不改用事后阈值 |
| GSE276045 的增殖/刚度结构不能分离 | 仅保留描述性方向，不能称为排除增殖混杂 |
| GWAS 版本、LD 或性状定义不可追溯 | 跳过 GWAS；不影响 S1/S2/S3 核心结果汇总 |
| 某一分支失败或结果阴性 | 作为边界证据写入，不追加无预设的新数据集 |
| 多重检验或模块覆盖不足 | 报全量结果和 `PARTIAL_COVERAGE`；不得放宽阈值 |

本轮执行完成条件是：S1–S3 每个分支都有“支持、未支持或不可估计”的明确裁决和文件化证据，而不是必须得到阳性结果。

## 7. 时间线与交付物

| 阶段 | 建议时间 | 交付物 |
|---|---:|---|
| 17A 方案锁定 | 第 0 周 | 本方案、claim contract、execution checklist |
| 17B–17C S1 | 第 1–2 周 | 来源/字段审计、样本内置换、共表达层级结果 |
| 17D S2 | 第 2–3 周 | GSE338388 provenance audit、双因素模块/候选统计 |
| 17E S3 | 第 3–5 周 | GSE123100 剂量–反应、GSE276045 竞争复核、可选 GWAS |
| 17F 综合 | 第 6 周 | evidence synthesis、claim revision map、论文修订动作 |

所有新增数据仍按项目目录保存：

`./data/metadata/independent_sources/`

`./data/processed/independent_sources/`

结果统一保存至：

`./results/14_second_round_computational_strengthening/`

## 8. 本轮最终决策

**Decision: PROCEED_WITH_LOCKED_SCOPE_AND_STRICT_BOUNDARIES**

执行方向已经确定：先 S1，后 S2，再做以 GSE123100 为主的 S3；GWAS 作为可选正交分支；不再扩大数据源范围，不再以计算阳性自动升级证据等级。

下一步是运行 **Step 17B：S1 来源与可分析性审计**。对应脚本为：

`./scripts/17B_audit_S1_source_feasibility.R`

该脚本仅读取已有审计结果和配置文件，不下载新数据、不初始化 Python、不进行假设检验。执行结果回来后再决定是否进入 S1 正式共表达分析。

## 9. Material Passport

- 输入：第二轮补强方案原文、当前证据状态、Step 12B/12C、Step 15F。
- 关键修正：将 GSE123100 锁定为 S3 正式剂量–反应主力，将 GSE276045/46 定位为 WI-38 刚度×增殖复核；将 S1 结论从 cell-intrinsic 判定降为共表达层级证据；将 S2 降为调控轴交叉验证。
- 输出：正式锁定版执行方案。
- 未执行内容：本文件未新增表达分析、未下载大文件、未改变 frozen panel 或当前证据等级。
- 下一执行节点：运行 `scripts/17B_audit_S1_source_feasibility.R`，完成 S1 provenance/可分析性审计。
