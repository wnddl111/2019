if (! requireNamespace("BiocManager"), quiet=TRUE)
  install.packages("BiocManager")
BiocManager::install("DESeq2") # :: biocManager�ȿ� �ִ� object�� �Լ��� ���ڴٴ� �� 

library(DESeq2)
library(ggplot2)

countdata_raw <- read.table("all.txt",header=TRUE, row.names=1) #header = TRUE��� �ϸ� HEADER�� ���̸��̾ߤ�~
# row.names =1 �̷��� 1���� row�̸����� �����ϰڴٴ� ��, r�� 1������ index��
head(countdata_raw) # r�� ǥ ���� �� �� �������� �׷��� �� ������������ Ȯ���Ϸ��� �տ��� 6������ �������� head�� Ȯ�� 
colnames(countdata_raw)
countdata <- countdata_raw[ ,6:ncol(countdata_raw)]
# countdata_raw���� ���� 5���� ������ ������ count���� �ƴϴ� �׷��� �̷��� ������
head(countdata)
colnames(countdata) # ���� colnames(���)=c("new���̸�1","~���̸�2",...) �̷��� ��/ �׳� �������ϸ� �����̸��� ������ 

colnames(countdata)=c("D217-1","D217-2",'Dj-1','Dj-2','Dn-1','Dn-2')
head(countdata)

######################################### ��������� ���� �ҷ��� column�� �ٲٱ� 


countdata_ko <- countdata[,c(1,2,5,6)] # ���ϰ� ���� �� ������ ����
head(countdata_ko)

count_not_ko <- countdata_ko # �ϳ��� ���Ժ��� ���ҰŰ�
count_norm_ko <- countdata_ko # �ϳ��� ���Ժ��� �ҰŴ� 


# transform countdata to matrix
count_norm_ko<-as.matrix(count_norm_ko)
count_not_ko <- as.matrix(count_norm_ko)

# 0 to NA
count_not_ko[count_not_ko==0] <- NA  # �̰� ���ϸ� ��հ��� 0���� �������ϱ� �װ� ����/ ��� �� ��ü ���� �������� ��

######## normalization �� �Ϸ���

# 1. coldata�� �������� coldata�� �� ���ܵ��� ������� ���� �� ���� 

condition <-factor(c("KO", "KO", "WT", "WT")) # ������ 2���� �����ϱ� ���� condition�̶�� ���� ���� 
coldata <- data.frame(row.names=colnames(count_norm_ko),condition)
head(count_norm_ko)
dds <- DESeqDataSetFromMatrix(count_norm_ko,coldata, ~ condition )
head(dds)
# deseqdataset~�̰� ���ĸ� wt1 wt2 ko1 ko2��� ���ܵ��� �־� �츮�� �̰� �����ں� �� �� ���� �󸶳� ���ǹ��� ���̰�
# �ִ��� ���ϰ� �;� �׷��� �� ������ wt�� ko�� ������ �װ� coldata�� ����ž�
# ���⼭ coldata�� condition �� wo,ko�� ���� �����ִ� count_norm_ko�� ������ �񱳸� �ϴ� � �⺻ �ڷḦ ����� ��  
dds<-dds[rowSums(counts(dds))>10,]
#���Ѵ� ���� 
DESeq_normed_ko <- DESeq(dds)
# �̰Ŵ� �׷��� dds�� ������estimating size factors,estimating dispersions ,gene-wise dispersion estimates
# mean-dispersion relationship, final dispersion estimates,fitting model and testing
# �̰͵��� ����ؼ� �� �࿡ ������ �߰��Ѵٴ� �� 
Result_ko <- results(DESeq_normed_ko, contrast = c("condition","KO","WT"))
# RESULT��� �Լ��� ���ݱ��� ���� �͵��� WT�� KO�� (�����߿�) KO/WT WT�� ���� KO�� ����� ���� 

# save DE analysis result
write.csv(Result_ko, "OsGASD-fch2-ko.csv")


# Draw Boxplot
dnorm_count_ko <- counts(DESeq_normed_ko, normalized=TRUE)
# float�� �ִ� ǥ�� ������� 
typeof(dnorm_count_ko)
# confirm normalization
head(dnorm_count_ko)

# 0 to NA
dnorm_count_ko[dnorm_count_ko==0] <- NA

########## visualization in one screen
dev.off() #�ٸ� �׸� ���� ���� ������ ������� �� 
par(mfrow = c(1,2), mar=c(6,5,3,1))
#par�� �ȷ�Ʈ ���� �Ű� mfrow�� �࿭ mar�� ������ ���� �ɷ� �������� ��ġ ���� 
#1�� 2��?

boxplot(log(count_not_ko+1,2) ~ col(count_not_ko), names = colnames(count_not_ko) ,cex.axis = 0.7, main = "before normalization", ylab = "log2(read counts + 1)" ,xlab="",las=0)
boxplot(log(dnorm_count_ko+1,2) ~ col(dnorm_count_ko), names = colnames(dnorm_count_ko) ,cex.axis = 0.7,main = "After normalization", ylab = "log2(read counts + 1)" ,xlab="", las=0)
# log 2�� �ϴµ� +1 �ϴ� ������ 0�ΰ��� ��� ����� �ȳ����ϱ� 1�� ���ؼ� 0���� ����� �ذ�, main�� ���� ����, ylab�� y�� ����, CEX.AXIS ->Y�� ǥ�� ũ��
# LAS=2 �ϸ� X�� LABLE�� ���η� ���ڰ� ���η� ��µǵ��� �ϴ� �� 
#�������� ������?
######################################### OX vs WT
countdata_raw = read.table("all.txt",header = TRUE, row.names=1)
head(countdata_raw)
countdata <- countdata_raw[ ,6:ncol(countdata_raw)]
head(countdata)
colnames(countdata)=c("KO","OSGASD-M",'OSGASD-OX','OX','WT-1st','WT-2nd')
head(countdata)


countdata_ox = countdata[, 3:6]
head(countdata_ox)
count_norm_ox <- countdata_ox
count_not_ox <- countdata_ox


count_norm_ox= as.matrix(count_norm_ox)
count_not_ox=as.matrix(count_not_ox)

count_not_ox[count_not_ox == 0] <-NA

condition = factor(c("OX","OX","WT","WT"))
coldata = data.frame(row.names = colnames(count_norm_ox),condition)

#1/0ó�� ���Ѵ밡 �Ǵ� ���� ���� �ϱ� ���� �� - �ι�° �� 
#�������� ref�� 0�ε� ���� �ϳ��� ���Ծ� �ƴ� �� 0.1�� �� �̷��� ������ ���Ѵ�� �����ϱ�

dds = DESeqDataSetFromMatrix(count_norm_ox,coldata,~ condition) # ~ �� ~�� ����
dds<-dss[rowSums(counts(dds))>10,]
DESeq(dds) -> DESeq_normed_ox

Result_ox <- results(DESeq_normed_ox, contrast = c("condition","OX","WT"))

write.csv(Result_ox,"OsGASD-fch2-ox.csv")
dnorm_count_ox <- counts(DESeq_normed_ox, normalized=TRUE)
head(dnorm_count_ox)

dnorm_count_ox[dnorm_count_ox==0] <- NA

par(mfrow = c(1,2), mar=c(5,5,3,1))

colnames(count_not_ox)-> xname

boxplot(log(count_not_ox+1,2) ~ col(count_not_ox), names = xname ,cex.axis = 0.7, main = "before normalization", ylab = "log2(read counts + 1)" ,xlab="",las=2)
boxplot(log(dnorm_count_ox+1,2) ~ col(dnorm_count_ox), names = colnames(dnorm_count_ox) ,cex.axis = 0.7,main = "After normalization", ylab = "log2(read counts + 1)" , xlab="",las=2)



head(count_not_ox)
colnames(count_not_ox)
head(dnorm_count_ox)
