import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 读取元数据文件
meta_data = pd.read_csv('g:/BB/meta.csv')

# 数据预处理
# 移除p值或logratio为NA的行
meta_data = meta_data.dropna(subset=['p value', 'logratio'])

# 计算-log10(p value)和-log10(FDR)用于y轴
meta_data['-log10(p)'] = -np.log10(meta_data['p value'])
meta_data['-log10(FDR)'] = -np.log10(meta_data['FDR'])

# 设置显著性阈值
p_threshold = 0.05
fc_threshold = 1.0  # log fold change阈值

# 创建火山图
plt.figure(figsize=(10, 8))
sns.set_style("whitegrid")

# 首先绘制所有点为灰色
plt.scatter(
    meta_data['logratio'],
    meta_data['-log10(p)'],
    alpha=0.5,
    color='grey',
)

# 定义显著上调和下调的代谢物
up_regulated = meta_data[(meta_data['logratio'] > fc_threshold) & (meta_data['p value'] < p_threshold)]
down_regulated = meta_data[(meta_data['logratio'] < -fc_threshold) & (meta_data['p value'] < p_threshold)]

# 绘制上调代谢物为红色
plt.scatter(
    up_regulated['logratio'],
    up_regulated['-log10(p)'],
    alpha=0.8,
    color='red',
    label=f'Up-regulated ({len(up_regulated)})'
)

# 绘制下调代谢物为蓝色
plt.scatter(
    down_regulated['logratio'],
    down_regulated['-log10(p)'],
    alpha=0.8,
    color='blue',
    label=f'Down-regulated ({len(down_regulated)})'
)

# 添加参考线
plt.axhline(y=-np.log10(p_threshold), color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=fc_threshold, color='gray', linestyle='--', alpha=0.5)
plt.axvline(x=-fc_threshold, color='gray', linestyle='--', alpha=0.5)

# 标记最显著的代谢物
# 选择前10个最显著的差异代谢物(按p值排序)
top_metabolites = pd.concat([up_regulated, down_regulated]).sort_values('p value').head(10)

for idx, row in top_metabolites.iterrows():
    plt.annotate(
        row['Metabolite'],
        (row['logratio'], row['-log10(p)']),
        fontsize=9,
        xytext=(5 if row['logratio'] > 0 else -5, 5),
        textcoords='offset points',
        ha='left' if row['logratio'] > 0 else 'right'
    )

# 设置轴标签和标题
plt.xlabel('Log2 Fold Change', fontsize=14)
plt.ylabel('-log10(p-value)', fontsize=14)

# 添加图例
plt.legend(fontsize=12)

# 保存图像
plt.tight_layout()
plt.savefig('g:/BB/metabolite_volcano_plot.png', dpi=300)
plt.show()