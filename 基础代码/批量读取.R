🟩1. 清空环境与加载核心包

rm(list = ls())    # 清空当前 R 工作区中的所有变量和对象，防止旧数据干扰
library(Seurat)    # 加载单细胞分析的核心 R 包——Seurat（目前使用的是 Seurat V5 版本）

🟩2. 批量读取数据并创建 Seurat 对象（循环部分）

f = dir("01_data/")       # 获取 01_data/ 目录下所有的子文件夹名称（即样本名列表：sample1, sample2...）
scelist = list()          # 创建一个空的列表（List），用来统一存放每个样本创建好的 Seurat 对象

for(i in 1:length(f)){    # 开始循环，遍历每一个样本
  # 1. 动态拼接路径并读取 10X 矩阵数据（barcodes, features, matrix）
  pda <- Read10X(paste0("01_data/", f[[i]])) 
  
  # 2. 将读取到的表达矩阵转化为 Seurat 对象，并以样本名作为项目标签（project）
  scelist[[i]] <- CreateSeuratObject(counts = pda, project = f[[i]])
  
  # 3. 打印当前样本的维度（基因数 和 细胞数）
  print(dim(scelist[[i]]))
}

📍用途：避免手动写重复代码。通过 for 循环，自动去每个样本文件夹里读数据。
💁输出解读：代码下方蓝色的数字（如 ## [1] 33538 10218）就是 print(dim(...)) 的结果：
33538 代表该样本中检测到的基因数量。
10218 代表该样本中包含的原始细胞数量。
你一共有 6 行输出，说明一共有 6 个样本。

🟩3. 合并数据集（Merge 与 JoinLayers）
# 1. 将第一个样本的对象与后面所有样本（scelist[-1]）合并成一个大对象
sce.all = merge(scelist[[1]], scelist[-1])
# 2. 整合层（Seurat V5 特有操作）
sce.all = JoinLayers(sce.all)

💬merge() 的用途：把分散在列表里的 6 个独立样本拼成一个叫 sce.all 的大矩阵，方便后续统一进行质控（QC）和降维。
💬JoinLayers() 的用途：这是 Seurat V5 的关键新特性。在 V5 中，merge 之后各个样本的表达矩阵默认是分开存储的（称为不同的 Layers，便于后续做高级整合）。这里运行 JoinLayers() 是把这些分开的矩阵重新融合成一个连续的、巨大的稀疏矩阵，以便于传统下游分析或计算。


🟩4. 查看细胞元数据（Meta-data）
head(sce.all@meta.data)

用途：查看合并后数据集前 6 行的细胞属性信息。
输出字段解读：
orig.ident：细胞来自哪个样本（如 sample1）。
nCount_RNA：该细胞内检测到的总 UMI（分子）数量（表达量总和）。
nFeature_RNA：该细胞内检测到的不同基因的数量。
左侧长长的字符串（如 AAACCCACACAAATAG-1_1）是每个细胞独一无二的条形码（Barcode），后面的 _1 是合并时自动加上的样本标记。

🟩5. 统计与校验
table(sce.all$orig.ident)

用途：统计大对象中，每一个样本分别包含多少个细胞。
输出解读：sample1 有 10218 个细胞，sample4 最多，有 12733 个细胞。这些数字与上面循环里打印出来的完全吻合。

sum(table(Idents(sce.all)))
用途：计算当前整张图谱里所有样本的细胞总数。
输出解读：## [1] 62348。也就是说，你目前已经成功把 6 个样本、共计 62,348 个细胞的数据全部打包进 sce.all 对象中了。



