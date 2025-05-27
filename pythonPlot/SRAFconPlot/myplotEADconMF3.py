import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

## 该函数输出的是AF情况下男性和女性离子电导缩放值的箱线图
Date = pd.read_csv('E:/pycharmProjects/py_Gender/A9SRAFconPlot/merged_file3.csv', header=None)

fig, ax =plt.subplots(2, 7, constrained_layout=True, figsize=(10, 6))
plt.tick_params(labelsize=7)
flatui = [[0/255, 104/255, 181/255], [242/255, 130/255, 38/255]]


pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 0], data=Date, ax=ax[0, 0], linewidth=0.8, palette=flatui,legend=False, fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[0, 0].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[0, 0].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 1], data=Date, ax=ax[0, 1], linewidth=0.8, palette=flatui,legend=False, fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[0, 1].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[0, 1].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 2], data=Date, ax=ax[0, 2], linewidth=0.8, palette=flatui,legend=False, fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[0, 2].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[0, 2].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 3], data=Date, ax=ax[0, 3], linewidth=0.8, palette=flatui,legend=False, fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[0, 3].set_ylim(bottom=0, top=4)  # 设置y轴范围
# 设置轴刻度值
ax[0, 3].set_yticks([0, 1, 2, 3, 4])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 4], data=Date, ax=ax[0, 4], linewidth=0.8, palette=flatui,legend=False, fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[0, 4].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[0, 4].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 17], y=Date.iloc[:, 5], data=Date, ax=ax[0, 5], linewidth=0.8, palette=flatui,legend=False, fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[0, 5].set_ylim(bottom=0, top=2)  # 设置y轴范围
# 设置轴刻度值
ax[0, 5].set_yticks([0, 1, 2])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 6], data=Date, ax=ax[0, 6], linewidth=0.8, palette=flatui,legend=False,fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[0, 6].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[0, 6].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 17], y=Date.iloc[:, 7], data=Date, ax=ax[1, 0], linewidth=0.8, palette=flatui,legend=False,fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[1, 0].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[1, 0].set_yticks([0, 1, 2, 3,4])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 8], data=Date, ax=ax[1, 1], linewidth=0.8, palette=flatui,legend=False,fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[1, 1].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[1, 1].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 17], y=Date.iloc[:, 9], data=Date, ax=ax[1, 2], linewidth=0.8, palette=flatui,legend=False,fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[1, 2].set_ylim(bottom=0, top=4)  # 设置y轴范围
# 设置轴刻度值
ax[1, 2].set_yticks([0, 1, 2, 3, 4])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 10], data=Date, ax=ax[1, 3], linewidth=0.8, palette=flatui,legend=False,fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[1, 3].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[1, 3].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 11], data=Date, ax=ax[1, 4], linewidth=0.8, palette=flatui,legend=False,fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[1, 4].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[1, 4].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 12], data=Date, ax=ax[1, 5], linewidth=0.8, palette=flatui,legend=False,fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[1, 5].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[1, 5].set_yticks([0, 1, 2, 3])

pic = sns.boxplot(x=Date.iloc[:, 16], y=Date.iloc[:, 13], data=Date, ax=ax[1, 6], linewidth=0.8, palette=flatui,legend=False,fliersize=0, width=0.7, saturation=6)
# 设置坐标轴范围
ax[1, 6].set_ylim(bottom=0, top=3)  # 设置y轴范围
# 设置轴刻度值
ax[1, 6].set_yticks([0, 1, 2, 3])

# 显示图形
plt.show()