RCDscoreR 使用说明（v0.1）
1) 包的核心功能
对输入表达矩阵（bulk / scRNA pseudo-bulk / 蛋白强度等）按每种 RCD 计算：
Driver-score：机制促发证据分数
Suppressor-score：抑制/对抗证据分数
Marker-score：表型指示证据分数
Net tendency = Driver - Suppressor（衍生指标，谨慎解释）
Conflict flag：该 RCD 的 driver 与 suppressor 基因集存在较高重叠（结构冲突）
同时引入 housekeeping 背景参照，使基因贡献只在“显著表达”时计入。

2) 输入数据要求（必须满足，否则结果不可靠）
2.1 输入表达矩阵 X
类型：数值矩阵（matrix / data.frame 可转 matrix）
维度：genes × samples/cells
rownames(X)：基因符号（必须与 genesets 使用的符号一致；人类建议全大写）
colnames(X)：样本/细胞名称（可缺省，会自动命名 S1,S2,…）
2.2 数据尺度（强烈建议）
由于算法依赖“以 housekeeping 为背景的 z-score”，X 建议是近似连续尺度：
bulk：log2(TPM+1) / logCPM / vst / rlog
scRNA：Seurat @data（log-normalized），或先做 pseudo-bulk 后再输入
protein：log2(intensity+1) / log2(LFQ+1) 等
不建议直接用 raw counts（必须先 log1p 或标准化）
1) 用 quantile 作为显著表达阈值的定义（你选的 A）
对每个样本/细胞 𝑗，计算阈值：
𝑇_𝑗=𝑄_𝑞 (𝑥_(⋅,𝑗) )，默认 q=0.75
定义基因贡献：
𝑎_(𝑔,𝑗)=max⁡(0ⓜ,𝑥_(𝑔,𝑗)−𝑇_𝑗 )
然后再对 driver/suppressor/marker 基因集做均值（或加权均值）聚合。
2) 对每类基因集聚合成 RCD 分数
对每个 RCD 𝑟：
Driver-score：
𝑆_(𝑟,𝑗)^𝑑𝑟𝑣=1/(∑128▒𝑤_𝑔 ) ∑128_(𝑔∈𝐺_𝑟^𝑑𝑟𝑣)▒𝑤_𝑔 ⋅𝑎_(𝑔,𝑗)
Suppressor-score（建议先独立输出，不要直接当负号）：
𝑆_(𝑟,𝑗)^𝑠𝑢𝑝=1/(∑128▒𝑤_𝑔 ) ∑128_(𝑔∈𝐺_𝑟^𝑠𝑢𝑝)▒𝑤_𝑔 ⋅𝑎_(𝑔,𝑗)
Marker-score：
𝑆_(𝑟,𝑗)^𝑚𝑘𝑟=1/(∑128▒𝑤_𝑔 ) ∑128_(𝑔∈𝐺_𝑟^𝑚𝑘𝑟)▒𝑤_𝑔 ⋅𝑎_(𝑔,𝑗)
Net tendency（可选，谨慎）：
𝑆_(𝑟,𝑗)^𝑛𝑒𝑡=𝑆_(𝑟,𝑗)^𝑑𝑟𝑣−𝑆_(𝑟,𝑗)^𝑠𝑢𝑝
4) Conflict flag（两层：基因集结构冲突 + 样本内冲突）
结构冲突（RCD 固有）：
𝑐𝑜𝑛𝑓𝑙𝑖𝑐𝑡\_𝑓𝑟𝑎𝑐_𝑟=(∣𝐺_𝑟^𝑑𝑟𝑣∩𝐺_𝑟^𝑠𝑢𝑝∣)/(∣𝐺_𝑟^𝑑𝑟𝑣∪𝐺_𝑟^𝑠𝑢𝑝∣)
超过阈值（比如 0.05 或 0.1）记为 TRUE。
样本内冲突（表达层面）：同一样本中 driver 与 suppressor 都显著激活时标记，例如：
𝑆_(𝑟,𝑗)^𝑑𝑟𝑣与 𝑆_(𝑟,𝑗)^𝑠𝑢𝑝都 > 各自全体分布的 75% 分位数（或固定阈值）
<img width="992" height="986" alt="image" src="https://github.com/user-attachments/assets/d54aff44-bb55-4a30-89d7-20b3428cf249" />
