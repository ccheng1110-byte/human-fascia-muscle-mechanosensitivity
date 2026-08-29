# 本轮研究补强过程记录

## 0. 文档信息

- 项目：人类筋膜–肌肉机械敏感性研究
- 记录范围：从补强意见采纳到 Step 15E 结果重综合、论文补强和 Step 16A 方案锁定
- 记录日期：2026-08-26
- 当前总体证据等级：**CAUTION**
- 当前阶段：计算补强已完成；功能验证方案已锁定；尚无新的湿实验结果

## 1. 补强启动时的研究状态

在补强启动前，研究已经完成了 GSE173252 发现分析、PRJNA607098 配对样本级验证、GSE130973 外部三角验证和证据等级综合。主要问题不是“缺少更多候选基因”，而是以下证据缺口：

1. Integrin/focal-adhesion 信号的竞争鲁棒性未通过，无法与 ECM 重塑、TGF/纤维化、炎症、缺氧和细胞周期程序充分区分。
2. Actomyosin 候选基因在独立来源中存在不一致，不能把单个基因或全部候选基因包装成保守机制。
3. PIEZO2 在不同来源中方向不一致，不能作为必需的机械感受器结论。
4. 既有数据主要是转录组关联，缺少机械扰动和功能表型。
5. GSE300230 的年龄与细胞系完全混杂，且只有两个培养成纤维细胞背景，不能支持独立年龄效应或人群供体泛化。
6. 既有培养细胞矩阵不能估计细胞内在调控与细胞组成变化的相对贡献。

因此，补强目标被重新定义为：**收窄关键不确定性、增加一个有明确机械对照的扰动证据层，并把最重要的湿实验验证问题具体化；而不是无限扩展公共数据或提高语言承诺。**

## 2. 意见采纳与范围取舍

### 2.1 采纳的核心意见

| 建议 | 采纳方式 | 采纳理由 |
|---|---|---|
| 引入公共机械/化学扰动数据 | 选择 GSE300230 作为唯一主力补强源 | 同时包含机械张力、TGF-beta 和两个人类真皮成纤维细胞背景，最直接回应“关联是否具有扰动 plausibility” |
| 保留 frozen candidate panel | 继续使用原 15 基因，不替换、不新增候选 | 防止补强阶段出现事后挑选或阈值放宽 |
| 将候选基因与模块结果分开 | 分别报告候选支持、模块方向和功能解释 | 当前模块级 actomyosin 信号比单基因一致性更稳定 |
| 进行竞争程序审计 | 将 cell-cycle、ECM、TGF/fibrosis、inflammation 等保留为正式竞争解释 | 直接回应 Integrin specificity gate 失败的问题 |
| 增加 score sensitivity 和 negative controls | 执行 z-mean、z-median、rank-mean、centered-mean、leave-one-gene-out 和表达分箱匹配阴性对照 | 检验模块方向是否依赖某一种打分方式或基因组成 |
| 降低因果性表述 | 将“机制成立”改为“扰动 plausibility/生物学合理性” | 公共培养细胞矩阵不能替代供体重复功能实验 |
| 把湿实验目标具体化 | 形成 Step 16A 供体重复、机械表型和通路扰动方案 | 计算补强的最终价值是确定最有信息量的下一项实验 |

### 2.2 暂不纳入主线的建议

以下方向在前期方案中被核验或提出，但本轮没有纳入主线：

- 多个额外刚度、YAP/TAZ、ATAC-seq 或 ChIP-seq 数据集的全面扩展；
- GWAS/MAGMA/LDSC/SMR 遗传学正交分析；
- scVI/scArches 的 F7–F8 标签映射；
- Human Protein Atlas 蛋白层面补强；
- 统一原始 FASTQ 重新处理。

暂不纳入不是因为这些方向没有价值，而是因为它们不能直接解决当前最核心的限制，且会引入新的组织来源、处理层级、数据准入和解释边界。按照“精简补强”原则，本轮优先完成一个设计清楚、可追溯、可解释的机械扰动数据源和稳健性审计。

## 3. 执行过程

### Step 15A：GSE300230 来源与可行性审计

**执行内容：**

- 审计官方 GEO 元数据和补充文件目录；
- 识别 24 个 GSM、2 个可解析的主要生物学背景；
- 下载两个小型处理后计数矩阵；
- 检查机械条件、TGF-beta 条件、年龄信息和 frozen panel 覆盖率。

**结果：**

- 最佳 frozen candidate coverage：15/15；
- 两个培养人真皮成纤维细胞系：GM08401 和 GM09503；
- 设计包含 tension/relaxed 与 TGF-beta absent/present；
- 初始自动 gate：`HOLD`。

**初始 HOLD 的原因：**

1. 解析器把 `NoTension` 错误识别为含糊的 tension 字段；
2. 数据集只有两个细胞系，不满足原先要求至少四个生物学单位的可行性规则。

输出：

- [Step 15A 可行性审计](./results/12_computational_strengthening/15A_GSE300230_provenance_feasibility_audit/GSE300230_step15A_provenance_feasibility_decision_v1.md)

### Step 15A2：门控裁决与边界修正

**采纳的修正：**

- 确认 `no tension` 与 `tension` 为两个不同机械水平；
- 将 16 个样本解释为每个细胞系内完整的 2×2 mechanical-tension × TGF-beta 设计；
- 允许在严格边界下继续，但不把两个细胞系当作人群供体重复。

**裁决：**

- 原始 gate：`HOLD`；
- 修正 gate：`PARTIAL_PROCEED_WITH_STRICT_BOUNDARIES`；
- 年龄与细胞系完全混杂，禁止估计独立年龄效应；
- 证据上限仍为 `CAUTION`。

输出：

- [Step 15A2 门控裁决](./results/12_computational_strengthening/15A2_GSE300230_gate_adjudication/GSE300230_step15A2_gate_adjudication_v1.md)

### Step 15B：分析契约冻结

在查看条件特异性表达结果前，冻结以下分析规则：

- 主要矩阵：`GSE300230_raw_counts_tensionTGFb.csv.gz`；
- GM08401 与 GM09503 分开分析；
- batch 作为阻断因素；
- 主要对比：无外源 TGF-beta 条件下 tension versus relaxed；
- TGF-beta 存在条件和 mechanical-by-TGF-beta interaction 作为次要分析；
- frozen 15 基因 panel、模块、竞争程序和多重校正规则不变；
- 不进行群体水平年龄推断。

输出：

- [Step 15B 契约冻结](./results/12_computational_strengthening/15B_GSE300230_contract_freeze/GSE300230_step15B_contract_freeze_decision_v1.md)

### Step 15C：机械化学扰动分析

**主要结果：**

| 分析对象 | GM08401 | GM09503 | 综合解释 |
|---|---:|---:|---|
| Actomyosin/Rho 模块 | Up，q = 5.071 × 10⁻⁹ | Up，q = 1.197 × 10⁻⁵ | 两个细胞系方向一致，最强的新证据层 |
| Integrin/focal-adhesion 模块 | Up，q = 2.666 × 10⁻⁴ | Up，q = 0.115 | 方向一致但统计强度不均一 |
| Frozen candidates 在两系均为正 | — | — | 8/15 |
| 两系均达到 q ≤ 0.05 的候选基因 | — | — | 0/15 |

**关键竞争结果：**

- Cell-cycle 模块在两条细胞系中均升高，q 值分别为 0.00225 和 5.12 × 10⁻¹²；
- 因此 actomyosin 的结果支持机械扰动 plausibility，但不能与增殖相关反应完全区分；
- Integrin 仍不能被表述为机械敏感性特异通路。

**Step 15C gate：** `PARTIAL_DIRECTIONAL_PLAUSIBILITY`。总体证据等级保持 `CAUTION`。

输出：

- [Step 15C 决策](./results/12_computational_strengthening/15C_GSE300230_mechanochemical_perturbation/GSE300230_step15C_mechanochemical_perturbation_decision_v1.md)
- [Step 15C 模块统计](./results/12_computational_strengthening/15C_GSE300230_mechanochemical_perturbation/GSE300230_frozen_module_camera_statistics_v1.csv)
- [Step 15C 候选统计](./results/12_computational_strengthening/15C_GSE300230_mechanochemical_perturbation/GSE300230_frozen_candidate_statistics_v1.csv)

### Step 15D：评分敏感性与阴性对照

**执行内容：**

- 比较 z-mean、z-median、rank-mean 和 centered-mean 四种模块评分；
- 进行 leave-one-gene-out 分析；
- 进行表达分箱匹配的随机阴性基因集对照；
- 审计细胞内在效应与细胞组成效应是否可估计。

**结果：**

- Core module direction stable under primary score：`TRUE`；
- Actomyosin 方向在两条细胞系、各种评分方式和 leave-one-gene-out 中均保持；
- Integrin 方向也保持，但其跨背景强度较不均一；
- 阴性对照没有消除 actomyosin 评分效应；
- 细胞内在效应与组成效应：`NOT_ESTIMABLE`；
- Cell-cycle 仍是必须保留的替代解释。

**Step 15D 不改变证据等级。** 该步骤证明的是计算评分的稳健性，不是细胞内在机制或因果性。

输出：

- [Step 15D 决策](./results/12_computational_strengthening/15D_GSE300230_score_sensitivity_and_negative_controls/GSE300230_step15D_score_sensitivity_and_negative_control_decision_v1.md)

### Step 15E：证据重综合与主张修正

**综合结论：**

> GSE300230 为培养人真皮成纤维细胞中的机械张力相关 actomyosin 转录响应提供了扰动 plausibility，但细胞周期竞争、供体泛化、年龄混杂、细胞组成效应和筋膜特异性仍未解决。

因此，研究主假设从“普遍的单基因机械敏感性机制”修正为：

> 人类筋膜样成纤维细胞状态与一个依赖情境的 ECM/integrin–actomyosin 转录程序相关；机械张力可能在培养的人真皮成纤维细胞中诱导 actomyosin-associated response。Integrin 反应不够均一，PIEZO2 证据依赖来源，且该反应可能部分耦合于增殖相关转录。

**Step 15E 结论：**

- Strengthening package：`COMPLETE`；
- Interpretation class：`ACTOMYOSIN_PERTURBATION_PLAUSIBILITY_WITH_PROLIFERATION_COMPETITION`；
- Overall evidence grade：`CAUTION`；
- 不授权把 manuscript evidence grade 升级。

输出：

- [Step 15E 证据重综合](./results/12_computational_strengthening/15E_evidence_resynthesis/Step15E_evidence_resynthesis_v1.md)
- [Step 15E 主张修订映射](./results/12_computational_strengthening/15E_evidence_resynthesis/Step15E_claim_revision_map_v1.csv)
- [Step 15E 论文修改动作](./results/12_computational_strengthening/15E_evidence_resynthesis/Step15E_manuscript_revision_actions_v1.md)

## 4. 论文与图表输出

Step 15E 结果已被纳入论文补强版：

- 摘要加入 GSE300230 的机械扰动结果；
- Methods 加入分层分析、TGF-beta 对照、batch 阻断和年龄混杂边界；
- Results 加入 actomyosin、integrin、cell-cycle 和敏感性分析结果；
- Discussion 将 actomyosin 定位为 perturbation plausibility，并明确增殖竞争；
- Conclusion 和 Data/code availability 同步更新；
- 原有 Figure 2–4 保留并继续嵌入 Word 版本；
- 修复页眉版本号和分类语句排版问题。

输出：

- [补强版 Markdown 初稿](./results/11_manuscript_preparation/14A_full_manuscript_draft/14A_full_manuscript_draft_v2_strengthened.md)
- [补强版 Word 初稿](./results/11_manuscript_preparation/14A_full_manuscript_draft/14A_full_manuscript_draft_v4_strengthened_with_figures.docx)
- [论文补强 QA 记录](./results/11_manuscript_preparation/14D_strengthened_manuscript_review/14D_strengthened_manuscript_review_v1.md)

Word 版共 16 页，已完成渲染和逐页视觉检查：未发现图像裁切、标题溢出、页码错误或参考文献溢出。

## 5. 从“意见”到“结果”的证据变化

| 维度 | 补强前 | 补强后 | 仍然不能宣称 |
|---|---|---|---|
| Actomyosin | 主要来自观察性模块/候选分析，外部候选基因不一致 | 两个培养真皮成纤维细胞背景中出现方向一致的机械张力相关模块响应，且评分稳健 | 供体重复的功能因果机制 |
| Integrin | 方向性较强但 competition gate 失败 | 仍为方向性优先候选，并得到部分机械扰动支持 | 机械敏感性特异性 |
| PIEZO2 | 来源依赖 | 继续保留为次要情境分支 | 普遍机械感受器机制 |
| 竞争解释 | 主要由 PRJNA607098 竞争鲁棒性失败提示 | Cell-cycle 在 GSE300230 中被直接显示为共同响应，替代解释更明确 | 已完成竞争程序解耦 |
| 稳健性 | 主要为既有数据流程稳健性 | 新增评分敏感性、leave-one-gene-out 和匹配阴性对照 | 细胞内在效应与组成效应 |
| 因果性 | 缺少机械扰动层 | 增加了培养细胞层面的 perturbation plausibility | 筋膜特异性、疼痛相关因果性 |

## 6. 当前结论与风险控制

### 可以保留的结论

- 在当前分析的培养人真皮成纤维细胞中，机械张力相关 actomyosin 转录响应具有生物学 plausibility。
- 筋膜样成纤维细胞状态与多组分 ECM/integrin–cytoskeletal 程序的关联仍具有研究价值。
- Integrin/focal-adhesion 是值得优先进行功能扰动验证的候选分支。

### 不能写入强结论的内容

- fascia-specific mechanosensitivity mechanism；
- validated mechanosensor 或 universal PIEZO2 involvement；
- pain/fibromyalgia causal pathway；
- independent human age effect；
- cell-intrinsic rather than composition-associated mechanism；
- 仅凭 GSE300230 宣称供体级独立重复。

### 补强没有解决的问题

1. 仍只有两个培养细胞背景；
2. 年龄与细胞系完全混杂；
3. 细胞周期与 actomyosin 响应可能耦合；
4. Integrin competition robustness 仍未通过；
5. 没有牵引力、收缩力或其他正交功能表型；
6. 没有真实筋膜组织或疼痛表型层面的验证。

## 7. 补强后的下一步

已据此生成 Step 16A 功能验证方案：

- [Step16A 功能验证方案](./results/13_functional_validation_planning/16A_protocol_lock/Step16A_functional_validation_protocol_v1.md)
- [Step16A 终点矩阵](./results/13_functional_validation_planning/16A_protocol_lock/Step16A_endpoint_matrix_v1.csv)
- [Step16A 设计模板](./results/13_functional_validation_planning/16A_protocol_lock/Step16A_experimental_design_template_v1.csv)
- [Step16A 决策记录](./results/13_functional_validation_planning/16A_protocol_lock/Step16A_protocol_lock_decision_v1.md)

当前操作性状态为：**方案锁定，等待实验平台和供体材料确认；尚未产生功能验证结果。**

## 8. Material Passport

- 意见来源：中期评估、已精读文献解题映射、GSE273293 替代方案以及 `14A_干实验补强方案.md`。
- 主要输入：Step 12B 方向锁定、Step 12C 功能验证规格、Step 15A–15D 审计与结果、Step 15E 证据重综合。
- 本轮新增数据：GSE300230 官方处理后矩阵及其审计结果；没有下载 617 GiB 级别的完整原始数据，也没有下载全量 H5AD。
- 主要转换：来源审计 → 分析契约冻结 → 分层机械扰动分析 → 稳健性审计 → 证据重综合 → 论文主张修订 → 功能验证方案锁定。
- 证据完整性：补强结果已纳入，但总体等级仍为 `CAUTION`，所有未解决边界已保留。
