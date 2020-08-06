library(ggplot2)
library(RColorBrewer)


ana_GO <- read.table("gtd4-�������.txt", sep="\t", header = T)
#���� ���� �� ����
head(ana_GO)
#View(ana_GO)

ana_GO$Cluster <- factor(ana_GO$Cluster, levels=c('gtd4-up','gtd4-down'))
ana_GO$GO.Name <- factor(ana_GO$GO.Name, levels=rev(unique(ana_GO$GO.Name)))

b <- ggplot(ana_GO, aes(x = Cluster, y = GO.Name))
b <- b + geom_point(aes(color= Fold_enrichment,size = -log10(Hyper.p.value)))
b <- b + scale_colour_gradientn(colours = c('blue','purple','red'),limit = c(2,25),na.value = 'red')
#limit�� �Ѿ�� ���� na.value�� �ǰ� �̰� red�� ǥ���ϰڴٴ� �ǵ� ���� �̰� �ȵ��ư����� �ִ� �׷� �Ʒ��ɷ� �Ѵ�. ������ ����� ��

#b <- b + scale_colour_gradientn(colours = c('blue','purple','red'), values = c(2,10,25),breaks =c(2,10,25), rescaler = function(x,...) x, oob = identity, na.value = 'red')
b <- b + theme(axis.text.x = element_text(size = 10,angle = 45, hjust = 1, vjust = 1))
#�̰Ŵ� x�� y�࿡ ǥ�õǴ� ���� ��ġ ���� angle�� ������ hjust�ϸ� x �� �������� 1��ŭ ��������v�� ����������� 
#b <- b + facet_grid(Cluster~., scales = "free_y", space="free")
b
