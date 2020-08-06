if (! requireNamespace("BiocManager"), quiet=TRUE)
  install.packages("BiocManager")
BiocManager::install("DESeq2") # :: biocManager안에 있는 object인 함수를 쓰겠다는 말 

library(DESeq2)
library(ggplot2)

countdata_raw <- read.table("all.txt",header=TRUE, row.names=1) #header = TRUE라고 하면 HEADER가 열이름이야ㅑ~
# row.names =1 이러면 1열을 row이름으로 지정하겠다는 말, r은 1번부터 index함
head(countdata_raw) # r은 표 같은 걸 잘 못보여줌 그래서 잘 가져와졌는지 확인하려고 앞에서 6개정도 가져오는 head로 확인 
colnames(countdata_raw)
countdata <- countdata_raw[ ,6:ncol(countdata_raw)]
# countdata_raw에서 앞의 5개는 유전자 정보지 count값이 아니다 그래서 이렇게 수정함
head(countdata)
colnames(countdata) # 보통 colnames(행렬)=c("new열이름1","~열이름2",...) 이렇게 씀/ 그냥 지정안하면 현재이름들 보여줌 

colnames(countdata)=c("D217-1","D217-2",'Dj-1','Dj-2','Dn-1','Dn-2')
head(countdata)

######################################### 여기까지는 파일 불러서 column명 바꾸기 


countdata_ko <- countdata[,c(1,2,5,6)] # 비교하고 싶은 두 집단을 설정
head(countdata_ko)

count_not_ko <- countdata_ko # 하나는 정규분포 안할거고
count_norm_ko <- countdata_ko # 하나는 정규분포 할거다 


# transform countdata to matrix
count_norm_ko<-as.matrix(count_norm_ko)
count_not_ko <- as.matrix(count_norm_ko)

# 0 to NA
count_not_ko[count_not_ko==0] <- NA  # 이걸 안하면 평균값이 0으로 떨어지니까 그거 방지/ 행렬 내 전체 값을 기준으로 함

######## normalization 을 하려면

# 1. coldata를 만들어야해 coldata는 이 집단들을 어떤식으로 구분 할 건지 

condition <-factor(c("KO", "KO", "WT", "WT")) # 집단을 2개로 구분하기 위해 condition이라는 변수 만듬 
coldata <- data.frame(row.names=colnames(count_norm_ko),condition)
head(count_norm_ko)
dds <- DESeqDataSetFromMatrix(count_norm_ko,coldata, ~ condition )
head(dds)
# deseqdataset~이게 뭐냐면 wt1 wt2 ko1 ko2라는 집단들이 있어 우리는 이게 유전자별 즉 행 별로 얼마나 유의미한 차이가
# 있는지 비교하고 싶어 그래서 이 집단을 wt와 ko로 나눴고 그게 coldata를 만든거야
# 여기서 coldata를 condition 즉 wo,ko에 따라서 원래있던 count_norm_ko를 나눠서 비교를 하는 어떤 기본 자료를 만드는 것  
dds<-dds[rowSums(counts(dds))>10,]
#무한대 방지 
DESeq_normed_ko <- DESeq(dds)
# 이거는 그러한 dds를 가지고estimating size factors,estimating dispersions ,gene-wise dispersion estimates
# mean-dispersion relationship, final dispersion estimates,fitting model and testing
# 이것들을 계산해서 그 행에 정보를 추가한다는 것 
Result_ko <- results(DESeq_normed_ko, contrast = c("condition","KO","WT"))
# RESULT라는 함수는 지금까지 만든 것들을 WT대 KO로 (순서중요) KO/WT WT에 비해 KO가 어떤지를 말함 

# save DE analysis result
write.csv(Result_ko, "OsGASD-fch2-ko.csv")


# Draw Boxplot
dnorm_count_ko <- counts(DESeq_normed_ko, normalized=TRUE)
# float이 있는 표를 만들었군 
typeof(dnorm_count_ko)
# confirm normalization
head(dnorm_count_ko)

# 0 to NA
dnorm_count_ko[dnorm_count_ko==0] <- NA

########## visualization in one screen
dev.off() #다른 그림 파일 열려 있으면 닫으라는 말 
par(mfrow = c(1,2), mar=c(6,5,3,1))
#par은 팔레트 같은 거고 mfrow는 행열 mar는 마방진 같은 걸로 동서남북 위치 지정 
#1행 2열?

boxplot(log(count_not_ko+1,2) ~ col(count_not_ko), names = colnames(count_not_ko) ,cex.axis = 0.7, main = "before normalization", ylab = "log2(read counts + 1)" ,xlab="",las=0)
boxplot(log(dnorm_count_ko+1,2) ~ col(dnorm_count_ko), names = colnames(dnorm_count_ko) ,cex.axis = 0.7,main = "After normalization", ylab = "log2(read counts + 1)" ,xlab="", las=0)
# log 2로 하는데 +1 하는 이유는 0인값인 경우 결과가 안나오니까 1을 더해서 0으로 만들어 준것, main은 위에 제목, ylab은 y축 제목, CEX.AXIS ->Y축 표시 크기
# LAS=2 하면 X축 LABLE을 세로로 글자가 세로로 출력되도록 하는 것 
#나머지는 무슨말?
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

#1/0처럼 무한대가 되는 값을 제거 하기 위한 것 - 두번째 줄 
#발현량이 ref는 0인데 나는 하나가 나왔어 아님 뭐 0.1이 떠 이러면 갯수가 무한대로 잡히니까

dds = DESeqDataSetFromMatrix(count_norm_ox,coldata,~ condition) # ~ 은 ~를 따라서
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

