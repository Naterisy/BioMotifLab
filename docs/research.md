# BioMotifLab 研究论证与方法决策记录

检索日期：2026-11-16。以下为文献事实与项目设计，不代表实验结果。实验结果以实际运行的配置、原始日志与汇总表为准。

## 具体生物计算问题

输入一组可能受到同一转录因子调控的未比对 DNA 序列，发现一个固定宽度、无间隙且可解释的序列模式，并评估其在留出序列中的匹配是否超过指定随机背景。目标是说明序列间 GC 组成异质性如何影响短基序发现，而非宣称只凭序列便可证明转录因子结合或基因调控。CRP 的短启动子示例提供真实序列案例；合成数据提供真实序列不能轻易获得的已知位点标签。

EM 的潜在位点模型、Gibbs 采样、PWM/PSSM、背景校正与 Sequence Logo 均为已有研究。Lawrence–Reilly（1990）给出早期 EM 位点发现；Bailey–Elkan（1994）发展混合模型的多基序发现；Lawrence 等（1993）提供 Gibbs 搜索框架。项目应把这些作为理论与基线来源，而不能将独立复现称为原创理论。[1–3]

## 三个候选方案比较

评分为本项目在本科课程资源限制下的工程决策，5分为更合适；“成本”“难度”分数越高代表越经济/容易落地，不是客观领域排名。

| 候选方案 | 生物学价值 | 可验证增量创新 | 计算成本适宜性 | 实现可控性 | 实验可验证性 | 总分 |
|---|---:|---:|---:|---:|---:|---:|
| A：组成感知的收缩背景与序列覆盖种子筛选，统一双链 EM | 4 | 3 | 5 | 5 | 5 | 22 |
| B：多宽度候选发现与重采样稳定性排序 | 4 | 3 | 3 | 3 | 4 | 17 |
| C：CNN 结合位点预测与梯度/突变归因 | 5 | 3 | 2 | 2 | 3 | 15 |

A 的可检验假设是：全局碱基频率不能表达不同启动子之间的组成差异；逐序列背景向全局收缩，可以在估计方差与局部适配之间折中。对种子按“出现于多少条序列”计数，能降低单条序列重复片段对初始化的支配。两者是针对本任务提出并独立实现的方法设计，必须分别消融；若只有全流程胜出，不能把收益归因于任意单个组件。

B 将固定宽度问题扩展为宽度选择，能避免用户事先知道长度；但候选空间越大，选择偏差越重。若每个宽度均需重采样且重新发现，成本会增加；课程项目可把它作为后续研究而非粗糙地混用不同长度的原始似然。

C 参照 DeepBind 的序列特异性学习思路，适合大规模实验数据，但本次可用 CRP 样本太少。以几条启动子训练 CNN 容易过拟合；直接调用预训练模型又不构成自主算法研究。因此本阶段不采用该路线。[15]

**选择 A。** 优先得到可核验的适用条件、失败案例及统计流程。最终实现采用零阶、链对称的逐序列收缩背景；一阶或更高阶模型仅作为后续扩展，不能写成已经实现。已有 BioProspector 与 MotifSampler 已研究高阶背景，Kim 等又讨论上下文，所以“自适应/Markov背景”这一大概念不具备本项目可主张的首创性。[4–6]

## 算法贡献边界与对照

实施模型为双链 OOPS：每条序列恰有一个候选位点；全部有效窗口与其反向互补共同参与潜在位点后验归一化。N 被保留，含 N 的窗口跳过，绝不能把未知碱基强行填成随机碱基。坐标统一为原始正链零基、右开区间，负链匹配单独给出 strand。

设序列 i 的链对称计数为 c_ib，全局平滑频率为 g_b，收缩强度 τ，使用 b_ib=(c_ib+τg_b+1/2)/(Σ_b c_ib+τ+2)。其模型假设仍为碱基独立，故不能校正二核苷酸依赖或局部短片段重复；这应成为失败条件而非隐藏因素。τ=0 是无跨序列收缩的局部估计消融，τ 很大趋于全局背景。

候选词以正反互补等价类聚合，并按每条序列至多一次的覆盖计数。以 Poisson 近似计算背景下至少出现一次的概率，得到覆盖偏离量，选择有限个初始化候选。这里的近似只用于启发式候选排序；因为重叠词并不独立，它不是正式 p 值。所有 EM 基线保持宽度、起点上限、停止准则和 PWM 伪计数一致。比较应包含：均匀背景 EM、全局背景 EM、收缩局部背景 EM、单 k-mer 种子对照。k-mer 对照未执行 EM，说明其是简单参考而非与成熟 MEME 等价的强基线。

消融至少包含固定全局背景、τ=0、随机种子替换覆盖富集种子。对同一合成数据种子配对比较，报告均值与标准差及各次明细。以 CPU 时间与运行峰值内存说明预算；不要通过给本方法更多重启来制造收益。若加入 MEME/STREME 实验，必须真实安装、记录版本和命令，不能将自己的简化 EM 命名为“MEME实测”。

## 统计有效性与评测协议

发现训练集与检验集按整条序列分开，避免窗口随机分割使几乎相同的序列泄漏到两侧。留出序列的统计量为双链所有有效窗口中的最大对数似然比，背景复制应使用相同长度、有效窗口/缺失模式和同一评分流程。在选择最佳窗口后仍用单个窗口零分布会高估显著性。

若采用逐序列单碱基打乱，检验的是在给定碱基计数、可交换排列假设下的富集；它保留 GC 组成但破坏二核苷酸与长程相关。若使用估计概率直接采样，则是参数自助零分布，p 值有背景估计误差。二者不能不加区别地称为精确检验。用 p=(1+超越次数)/(B+1)，报告 B 和最小可分辨 p，不能给出零 p。[13]

固定训练好的 PWM 后在独立留出序列检验，能够避开对同一数据“先挑基序、再给显著性”的一部分选择偏差；训练集分数只能用于描述。一次检验多个留出序列时输出 BH 调整值并说明其独立/相关条件，不把 q≤0.05 解读为每个位点有95%概率真实结合。[14]

合成实验同时保存植入位点、正负链、随机种子和生成分布；定位 Precision/Recall/F1 应明确是核苷酸重叠、窗口IoU还是位点容差，而不是只比较共识字符串。允许无基序序列时，OOPS 强制定位的模型错误必须计为限制；不能把每条序列都有植入位点的实验结果推广为有无结合位点的分类效果。真实 CRP 案例无导入的核实定位标签，故只报告模型恢复、留出打乱对照和稳定性，不能报告生物学定位 F1。真实数据注释不完备与合成背景偏差是经典基准已强调的问题。[12]

## 数据、工具与许可核实

MEME Suite 是成熟的基序工具集合，覆盖发现、已知基序扫描及比较。[7–9] 它不是本项目需要“替代”的落后工具。课程项目的价值是透明公式、可复现实验与局部组成适配的可解释研究；不应无证据宣传精度或速度超过 MEME/STREME。官方文档 <https://meme-suite.org/meme/doc/meme.html> 提供 CRP 的18条105nt真实示例及固定宽度20的命令。官方下载 URL、字节哈希、大小和使用范围见 `data/real/manifest.json`。官方许可证 <https://meme-suite.org/meme/doc/copyright.html> 有教育/科研/非营利条件，不应把其源代码或数据整体重新授权为MIT。因此只公开数据下载器、校验与必要汇总，不公开其原始受约束文件。

JASPAR 2026 官方页面 <https://jaspar.elixir.no/> 明确列出该版本论文与 CC BY 4.0 授权。CORE 是实验支持的真核转录因子 PFM 集合；UNVALIDATED 的证据类型不同，应分开。固定 MA0139.1 CTCF profile 提供可复现的外部模型示例，API记录原始实验引用 PMID 17512414。PFM不是原始序列，也不与细菌CRP对应。已下载该矩阵并保留署名；从它采样生成的序列必须标为合成数据。[10]

Sequence Logo 使用每列熵与字母频率表达模式，其基本思想已有文献。[11] 标准基于均匀参考的2−H信息量图与非均匀背景下的对数比不同，应在图注说明，不把图上高度误解为检验显著性。

## 已核验文献

[1] Lawrence C E, Reilly A A. An expectation maximization (EM) algorithm for the identification and characterization of common sites in unaligned biopolymer sequences[J]. Proteins, 1990, 7(1):41–51. DOI: [10.1002/prot.340070105](https://doi.org/10.1002/prot.340070105). [核验来源](https://pubmed.ncbi.nlm.nih.gov/2184437/)。

[2] Bailey T L, Elkan C. Fitting a mixture model by expectation maximization to discover motifs in biopolymers[C]. Proceedings of the Second International Conference on Intelligent Systems for Molecular Biology, 1994, 28–36. 未填DOI（避免猜测）。 [核验来源](https://cdn.aaai.org/ISMB/1994/ISMB94-004.pdf)。

[3] Lawrence C E, Altschul S F, Boguski M S, Liu J S, Neuwald A F, Wootton J C. Detecting subtle sequence signals: a Gibbs sampling strategy for multiple alignment[J]. Science, 1993, 262(5131):208–214. DOI: [10.1126/science.8211139](https://doi.org/10.1126/science.8211139). [核验来源](https://www.cs.cmu.edu/~durand/03-711/Readings/gibbsSampler.pdf)。

[4] Liu X, Brutlag D L, Liu J S. BioProspector: discovering conserved DNA motifs in upstream regulatory regions of co-expressed genes[C]. Pacific Symposium on Biocomputing, 2001, 6:127–138. 未填DOI（避免猜测）。 [核验来源](https://psb.stanford.edu/psb-online/proceedings/psb01/liu.pdf)。

[5] Thijs G, Lescot M, Marchal K, Rombauts S, De Moor B, Rouzé P, Moreau Y. A higher-order background model improves the detection of promoter regulatory elements by Gibbs sampling[J]. Bioinformatics, 2001, 17(12):1113–1122. DOI: [10.1093/bioinformatics/17.12.1113](https://doi.org/10.1093/bioinformatics/17.12.1113). [核验来源](https://academic.oup.com/bioinformatics/article/17/12/1113/225238)。

[6] Kim N K, Tharakaraman K, Spouge J L. Adding sequence context to a Markov background model improves the identification of regulatory elements[J]. Bioinformatics, 2006, 22(23):2870–2875. DOI: [10.1093/bioinformatics/btl528](https://doi.org/10.1093/bioinformatics/btl528). [核验来源](https://academic.oup.com/bioinformatics/article/22/23/2870/279908)。

[7] Bailey T L, Johnson J, Grant C E, Noble W S. The MEME Suite[J]. Nucleic Acids Research, 2015, 43(W1):W39–W49. DOI: [10.1093/nar/gkv416](https://doi.org/10.1093/nar/gkv416). [核验来源](https://pmc.ncbi.nlm.nih.gov/articles/PMC4489269/)。

[8] Grant C E, Bailey T L, Noble W S. FIMO: scanning for occurrences of a given motif[J]. Bioinformatics, 2011, 27(7):1017–1018. DOI: [10.1093/bioinformatics/btr064](https://doi.org/10.1093/bioinformatics/btr064). [核验来源](https://pmc.ncbi.nlm.nih.gov/articles/PMC3065696/)。

[9] Bailey T L. STREME: accurate and versatile sequence motif discovery[J]. Bioinformatics, 2021, 37(18):2834–2840. DOI: [10.1093/bioinformatics/btab203](https://doi.org/10.1093/bioinformatics/btab203). [核验来源](https://academic.oup.com/bioinformatics/article/37/18/2834/6184861)。

[10] Baydar Ovek D, Rauluseviciute I, Aronsen D R, et al.. JASPAR 2026: expansion of transcription factor binding profiles and integration of deep learning models[J]. Nucleic Acids Research, 2026, 54(D1):D184–D193. DOI: [10.1093/nar/gkaf1209](https://doi.org/10.1093/nar/gkaf1209). [核验来源](https://jaspar.elixir.no/)。

[11] Schneider T D, Stephens R M. Sequence logos: a new way to display consensus sequences[J]. Nucleic Acids Research, 1990, 18(20):6097–6100. DOI: [10.1093/nar/18.20.6097](https://doi.org/10.1093/nar/18.20.6097). [核验来源](https://academic.oup.com/nar/article/18/20/6097/1141316)。

[12] Tompa M, Li N, Bailey T L, et al.. Assessing computational tools for the discovery of transcription factor binding sites[J]. Nature Biotechnology, 2005, 23:137–144. DOI: [10.1038/nbt1053](https://doi.org/10.1038/nbt1053). [核验来源](https://www.nature.com/articles/nbt1053)。

[13] Phipson B, Smyth G K. Permutation P-values should never be zero: calculating exact P-values when permutations are randomly drawn[J]. Statistical Applications in Genetics and Molecular Biology, 2010, 9(1):Article39. DOI: [10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585). [核验来源](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf)。

[14] Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing[J]. Journal of the Royal Statistical Society: Series B (Methodological), 1995, 57(1):289–300. DOI: [10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x). [核验来源](https://onlinelibrary.wiley.com/doi/full/10.1111/j.2517-6161.1995.tb02031.x)。

[15] Alipanahi B, Delong A, Weirauch M T, Frey B J. Predicting the sequence specificities of DNA- and RNA-binding proteins by deep learning[J]. Nature Biotechnology, 2015, 33:831–838. DOI: [10.1038/nbt.3300](https://doi.org/10.1038/nbt.3300). [核验来源](https://www.nature.com/articles/nbt.3300)。

