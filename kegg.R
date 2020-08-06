# bioconductor source 

# source("https://bioconductor.org/biocLite.R")

# install packages -> at first time 
#biocLite("clusterProfiler")
#biocLite("biomaRt")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#install.packages("stringi", type="source") #install�� �ѹ���! 
#install.packages("MASS")
#install.packages("devtools")
#devtools install


#BiocManager::install("https://guangchuangyu.github.io/software/clusterProfiler")
#devtools::install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))
#������ �������� ������ 
# attach packages

#clusterprifiler install
library("clusterProfiler")
library("ggplot2")

# set working directory
#setwd("C:\\Users\\�輮��\\Desktop")

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
#RAP_id�� cluster�� ���� �׸���. data�� ana�� �־��ִ¿� �� �Լ��� dataframe���� �Ȱ� ���ư��� �տ� ������ ������ϴ°�
# fun�� �Լ� �̸��̰� kegg enrich�� �ҰŶ� �̰� �־���, �׸��� organism�� �������� ������ �ΰ����� �ؾ����� ���ĸ��ؾ�����
#�𸣴ϱ� dosa�̰� ���� ���õȰǰ��� �̰� �ؾ���, pvalue�� 0.05���� �Ǵ� �� �ϰ� �ʹٴ� �Ҹ�
# Plotting Enrichment result
# because clusterProfiler used ggplot2 for its plot, we can add other graphics via ggplot2

xx@compareClusterResult$Cluster
# �̷��� ���� xx�� class�� �� �ȿ� ������ @����ؾ��ϰ� �� �ȿ� ���� �׾��� ���� ���� �ſ� ���� �Ϸ��� $�� ����.

xx@compareClusterResult$Cluster <- factor(xx@compareClusterResult$Cluster, levels=c('Up','Down'))
#���⼭ factor�� aa,bb�� ���� �� a,b�� factor���ϴµ� ���� �̷��� factor(��¼��, ���ñ�)�ϸ� ��¼���� ���ñ� ������ ���ڴٴ� ��
#�̰� ���� ���ϸ� ���ĺ� ������ �����ؼ� �ٲ��ذ�

kegg <- dotplot(xx) # normal dot plot, x tics overlapped
kegg <- kegg + facet_grid(.~Cluster, scales = "free_x", space="free") 
#dotplot�� �׸��׸��� �Լ��� ggplotó�� +�� �ɼ� �߰��� ��. facet_grid�� x,y ���� �̸�? �̷��Ŷ� ���õȰ� ������ 
#ù��° �Ķ���ʹ� yx������ ǥ���ϴµ� .�� y�ΰŰ� x�� ~cluster�� y�� �� �ǵ帮�� x�� �̰ſ� ���� �����ϰڴٴ� ��
# scales�� x�ǵ�����ϱ� x�� free�� �ƴ϶� free_x �̷��� �Ѱ�
# ��ǻ� y�� �� �ǵ���� space�� �� �ʿ� ����
kegg # plotting