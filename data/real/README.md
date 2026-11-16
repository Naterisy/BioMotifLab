# 真实数据与使用边界

检索与下载日期：2026-11-16。文件级来源、SHA-256、用途与许可见 `manifest.json`。

## MEME CRP 示例

官方下载地址：<https://meme-suite.org/meme-software/meme_example_inputs/crp0.s>；[官方命令示例](https://meme-suite.org/meme/doc/meme.html)将该文件用于 DNA 基序发现及双链搜索。下载的原文件含 18 条序列，每条 105 nt，共 1890 nt。历史 CRP 案例参见 Lawrence 与 Reilly（1990），DOI [10.1002/prot.340070105](https://doi.org/10.1002/prot.340070105)。

这是来源可核实的真实序列案例，适用于恢复基序、观察序列组成和留出对照检验；规模小、选择偏倚明显，不能代表全基因组性能。原 FASTA 标题后的数字原样保留，但未核实其坐标语义，不作为定位真值。本文若提供精确定位准确率，应来自有植入坐标的合成数据，不能来自对这组数据的猜测标注。把 MEME 的预测坐标当作实验真值也不合适。

[MEME 许可](https://meme-suite.org/meme/doc/copyright.html)允许教育、科研和非营利用途，附有保留完整版权声明等条件；这不等于无条件的商业开源许可。项目 MIT 许可不覆盖该数据。仓库仅保留下载脚本、来源清单与汇总统计；原 FASTA 和下载的完整版权页被忽略，不上传。使用者应阅读下载的 `MEME-COPYRIGHT.html` 并遵守其条款。

## JASPAR CTCF PFM

MA0139.1 是 JASPAR CORE 中的人类 CTCF 位置频率矩阵，19 列，profile 来源标注 ChIP-seq，关联 PMID 17512414。使用固定版本便于复现，不声称其为最新 CTCF profile。JASPAR 网站明确采用 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)；本项目保留原始 API JSON、版本与引用。署名：JASPAR contributors；Baydar Ovek D, et al. JASPAR 2026，DOI [10.1093/nar/gkaf1209](https://doi.org/10.1093/nar/gkaf1209)。

PFM 是实验导出的位点统计模型，不是未经处理的基因组序列。此文件用于已知模型读取/可视化或合成采样示例；它不是 CRP 的对应基序，也不能单独满足真实 DNA 序列实验要求。JASPAR CORE 主要收录真核转录因子，不应因其知名度将其硬套到细菌 CRP。

## 下载与完整性

```bash
python scripts/download_data.py
# 若系统要求显式 CA 包，例如当前 macOS 环境：
python scripts/download_data.py --cafile /etc/ssl/cert.pem
```

脚本使用验证 TLS，检查 CRP 文件 SHA-256 与 JASPAR 矩阵内容 SHA-256；上游内容变动时终止，不静默替换基准。API 中非生物学元数据可能更新，因此矩阵语义哈希与原始响应哈希分别记录。文件写入采用先临时文件、后原子替换。
