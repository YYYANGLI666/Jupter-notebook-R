#为 10x Genomics 单细胞原生数据进行重构与清洗

#在 GEO 数据库（比如这里的 GSE231920）下载单细胞原始数据时，所有样本的 (barcodes.tsv.gz)、(features.tsv.gz) 和 (matrix.mtx.gz )通常会被全塞在同一个文件夹里，并加上样本前缀（如 GSMxxxx_sample1_barcodes.tsv.gz）。
#而 Seurat 的 Read10X() 函数要求每个样本必须有独立文件夹，且里面三个文件的名字必须精确为 barcodes.tsv.gz、features.tsv.gz 和 matrix.mtx.gz。
#这段代码正是为了解决这个问题。

#第一部分：获取文件列表与提取样本名
#untar("GSE231920_RAW.tar",exdir = "GSE231920_RAW")  #解压原始tar包（已被注释）
#unlink("GSE231920_RAW.tar")                         #删除原始tar包（已被注释）
library(stringr)

# 1. 获取解压后文件夹内所有文件的完整路径
fs = paste0("GSE231920_RAW/",dir("GSE231920_RAW/"))
fs

# 2. 提取所有唯一的样本名
samples = dir("GSE231920_RAW/") %>% str_split_i("_",2) %>% unique();samples

#1+2：对前面代码的一个总结：
#fs：包含了所有文件的路径列表。
#samples：核心在 str_split_i("_", 2)。由于文件名类似于 GSM707_sample1_barcodes.tsv.gz，代码以 _ 作为分隔符进行切分，取第 2 个片段（即 sample1）。最后用 unique() 去重，从而拿到了所有样本名列表。

#第二部分：为每个样本创建独立文件夹
#为每个样本创建单独的文件夹
lapply(samples, function(s){
  ns = paste0("01_data/",s)
  if(!file.exists(ns))dir.create(ns,recursive = T)
})

#*遍历刚刚提取出的样本名列表（samples）。
#*在当前工作目录下创建一个叫 01_data 的新文件夹，并在其内部为每个样本建立一个以它名字命名的子文件夹（例如 01_data/sample1/）。
#*recursive = T 确保如果父文件夹 01_data 不存在，会自动连同子文件夹一起创建。

#第三部分：分发文件到对应的文件夹
#每个样本的三个文件复制到单独的文件夹
lapply(fs, function(s){
  #s = fs[1]
  for(i in 1:length(samples)){
    #i = 1
    if(str_detect(s,samples[[i]])){
      file.copy(s,paste0("01_data/",samples[[i]]))
    }
  }
})

#这是一个嵌套循环：外层遍历所有原始文件（fs），内层遍历所有样本名（samples）。
#str_detect(s, samples[[i]])：检查当前文件的路径字符串中，是否包含某一个样本名。
#如果匹配成功，就用 file.copy 把该文件复制到刚刚建好的对应样本文件夹中。

#第四部分：文件重名（关键步骤）
#⚠️ 教程警告：“注意这段改名的代码不要反复运行”。因为一旦名字已经改掉了，再次运行就匹配不到原来的前缀模式，反而会把正确的文件名搞乱。
#文件名字修改
on = paste0("01_data/",dir("01_data/",recursive = T));on
nn = str_remove(on,"GSM\\d+_sample\\d_");nn
file.rename(on,nn)

#on (Old Names)：获取 01_data/ 文件夹下所有复制过来的文件的当前路径（此时它们还带着长长的类似 GSM123456_sample1_ 的前缀）。
#nn (New Names)：使用正则表达式 GSM\\d+_sample\\d_ 匹配并删除这些前缀。

#GSM\\d+ 匹配 GSM 加上一串数字。
#sample\\d_ 匹配样本编号前缀。
#擦除前缀后，文件名就只剩下 barcodes.tsv.gz、features.tsv.gz 和 matrix.mtx.gz 了。
#file.rename(on, nn)：正式在本地磁盘上执行重命名。
#运行后的最终效果：
#你的工作目录下会多出一个 01_data 文件夹，结构如下，此时直接喂给 Seurat 的 Read10X("01_data/sample1") 就能完美读取了：
rm(list = ls())

