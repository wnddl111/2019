# bioconductor source 

# source("https://bioconductor.org/biocLite.R")

# install packages -> at first time 
#biocLite("clusterProfiler")
#biocLite("biomaRt")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#install.packages("stringi", type="source") #install은 한번만! 
#install.packages("MASS")
#install.packages("devtools")
#devtools install


#BiocManager::install("https://guangchuangyu.github.io/software/clusterProfiler")
#devtools::install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))
#위에거 에러나서 무시함 
# attach packages

#clusterprifiler install
library("clusterProfiler")
library("ggplot2")

# set working directory
#setwd("C:\\Users\\김석희\\Desktop")

## read clustered information with RAP id ## 

#    Classification          RAP_id
#1 Leaf + Flag leaf Os01t0105900-01
#2 Leaf + Flag leaf Os01t0208700-01

ana_rap <- read.table("KEGG.txt", header = T, sep = "\t")
View(ana_rap)

# convert matrix to dataframe
ana <- as.data.frame(ana_rap)

# Enrich KEGG of clustered data
xx <- compareCluster(RAP_id ~ Cluster , data = ana, fun='enrichKEGG', organism="dosa", pvalueCutoff=0.05)
head(xx)
#RAP_id를 cluster에 따라서 그린다. data는 ana를 넣어주는에 이 함수는 dataframe으로 된게 돌아가서 앞에 과정을 해줘야하는것
# fun은 함수 이름이고 kegg enrich를 할거라서 이걸 넣어줌, 그리고 organism은 정해주지 않으면 인간으로 해야할지 초파리해야할지
#모르니까 dosa이게 벼랑 관련된건가봐 이거 해야함, pvalue는 0.05정도 되는 걸 하고 싶다는 소리
# Plotting Enrichment result
# because clusterProfiler used ggplot2 for its plot, we can add other graphics via ggplot2

xx@compareClusterResult$Cluster
# 이렇게 만든 xx는 class임 그 안에 들어가려면 @사용해야하고 더 안에 들어가서 그안의 변수 같은 거에 접근 하려면 $를 쓴다.

xx@compareClusterResult$Cluster <- factor(xx@compareClusterResult$Cluster, levels=c('Up','Down'))
#여기서 factor는 aa,bb가 있을 때 a,b를 factor라하는데 지금 이렇게 factor(어쩌고, 뭐시기)하면 어쩌고를 뭐시기 순으로 보겠다는 것
#이걸 지정 안하면 알파벳 순으로 정렬해서 바꿔준것

kegg <- dotplot(xx) # normal dot plot, x tics overlapped
kegg <- kegg + facet_grid(.~Cluster, scales = "free_x", space="free") 
#dotplot도 그림그리는 함순데 ggplot처럼 +로 옵션 추가가 됨. facet_grid는 x,y 축의 이름? 이런거랑 관련된것 같은데 
#첫번째 파라미터는 yx순으로 표기하는데 .이 y인거고 x가 ~cluster라서 y는 안 건드리고 x는 이거에 따라서 구분하겠다는 것
# scales는 x건드렸으니까 x를 free가 아니라 free_x 이렇게 한것
# 사실상 y는 안 건드려서 space는 할 필요 없다
kegg # plotting
