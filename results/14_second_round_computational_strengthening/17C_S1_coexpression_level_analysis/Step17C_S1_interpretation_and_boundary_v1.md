# Step 17C S1 共表达层级结果解释与边界

## 中间裁决

**Interpretation class: SOURCE_DEPENDENT_MIXED_COEXPRESSION_LEVEL_SIGNAL**  
**Evidence grade: CAUTION**

Step 17C 完成了预设候选基因对的单细胞共检测、样本内置换和样本/受试者层面汇总。该步骤支持的是共表达层级证据，不是细胞内在因果机制。

## 主要结果

| 来源 | 统计单位 | 预设基因对 | 描述性支持 | 解释 |
|---|---:|---:|---:|---|
| PRJNA607098 | 12 个重建样本 | 54 | 37/54 | 目标 F7 状态内存在共表达层级信号 |
| GSE130973 | 5 个受试者 | 54 | 0/54 | 未达到预设 pair-level 支持 |

PRJNA607098 的 37/54 支持由以下配对家族组成：

- integrin × actomyosin：18/24；
- integrin × mechanosensor：10/18；
- integrin × hippo：9/12。

PRJNA607098 的目标状态细胞数为 30,803，与此前 Step 08B3 的 F7 标签结果一致。该数字对应 `F7: Fascia-like myofibroblast` 标签，不应与另一个 `Myofibroblast` 标签的 1,714 个细胞混用。

## 方法学限制

1. PRJNA607098 没有可用的 library-size 字段，因此本次置换保留了每个基因的样本内检出边际，但没有完成基于 library size 的分层校正。
2. PRJNA607098 的 SRS 号是样本级单位，不能自动等同于 12 个独立供体。
3. GSE130973 只有 5 个受试者，且为独立皮肤成纤维细胞来源；其阴性结果不能解释为机制不存在，但说明 PRJNA 信号没有在该来源中稳定复现。
4. 单细胞共检测受 dropout、测序深度、检测率和目标细胞定义影响。

## 允许的表述

> The S1 analysis identified source-dependent, mixed co-expression-level evidence for selected integrin–mechanotransduction candidate pairs. This pattern narrows the uncertainty about co-detection but does not establish a cell-intrinsic mechanism, causal mechanosensitivity, donor-independent replication, or fascia specificity.

## 不允许的表述

- `cell-intrinsic mechanosensitivity was proven`；
- `integrin mechanosensing is conserved across sources`；
- `PRJNA607098 provides 12 independent donor replications`；
- `the negative GSE130973 result disproves the program`。

## 对总体方案的影响

- Step 17C 不改变总体证据等级，继续保持 **CAUTION**。
- 不触发候选基因替换、阈值放宽或跨研究合并。
- S1 的价值是将“细胞内在 vs 组成效应”从完全 `NOT_ESTIMABLE` 收窄为“来源依赖的共表达层级信号”，但没有解决细胞内在因果问题。
- 下一步进入 Step 17D：GSE338388 的官方来源和实验设计审计，严格检验其是否真实支持 TGFβ 暴露（±）× TEAD 抑制（±）的调控轴交叉验证。

## Material Passport

- 输入：Step 17C PRJNA607098 与 GSE130973 共表达结果、Step 08B3 F7 规模核对结果。
- 转换：来源分层、pair-level 支持汇总和证据边界解释。
- 新数据下载：无。
- 证据状态：S1 完成，结论限制为来源依赖的共表达层级证据。
