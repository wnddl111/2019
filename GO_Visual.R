library(ggplot2)
library(RColorBrewer)


ana_GO <- read.table("gtd4-닷블랏용.txt", sep="\t", header = T)
#엑셀 파일 못 연다
head(ana_GO)
#View(ana_GO)

ana_GO$Cluster <- factor(ana_GO$Cluster, levels=c('gtd4-up','gtd4-down'))
ana_GO$GO.Name <- factor(ana_GO$GO.Name, levels=rev(unique(ana_GO$GO.Name)))

b <- ggplot(ana_GO, aes(x = Cluster, y = GO.Name))
b <- b + geom_point(aes(color= Fold_enrichment,size = -log10(Hyper.p.value)))
b <- b + scale_colour_gradientn(colours = c('blue','purple','red'),limit = c(2,25),na.value = 'red')
#limit가 넘어가는 경우는 na.value가 되고 이걸 red로 표시하겠다는 건데 가끔 이게 안돌아갈때가 있다 그럼 아래걸로 한다. 이유는 사수도 몲

#b <- b + scale_colour_gradientn(colours = c('blue','purple','red'), values = c(2,10,25),breaks =c(2,10,25), rescaler = function(x,...) x, oob = identity, na.value = 'red')
b <- b + theme(axis.text.x = element_text(size = 10,angle = 45, hjust = 1, vjust = 1))
#이거는 x축 y축에 표시되는 값의 위치 조절 angle은 각도고 hjust하면 x 축 기준으로 1만큼 떨어져라v는 수평기준으로 
#b <- b + facet_grid(Cluster~., scales = "free_y", space="free")
b

