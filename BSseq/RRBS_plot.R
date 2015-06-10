setwd("/lustre/user/liclab/zhangc/proj/lingte/RRBS/STAT")

###############  cumulative depth 
d3 = read.table("D3.depth.sort_r.acu.sort", row.names=1)
d4 = read.table("D4.depth.sort_r.acu.sort",row.names=1)

#d3$V2 = d3$V2/49235722
#d4$V2 = d4$V2/49349983


D3 = d3[0:100,]
D4 = d4[0:100,]
mat = cbind(D3, D4)

library(reshape2)
library(ggplot2)
m <- melt(mat)
colnames(m) = c("depth", "Sample", "Cumulative_Ratio")
p <- ggplot(m, aes(x=depth, y=Cumulative_Ratio, col=Sample)) + geom_line(size=2.0)
p +  theme( text = element_text(size=30),
            axis.line = element_line(size = 1),
            panel.background = element_blank()) +
 labs(title="Cumulative Depth of C sites ", x="Depth of C sites", y="Covered Bases")


############## methylation level

D3_CG = read.table("D3.methyl.CG")
D4_CG = read.table("D4.methyl.CG",row.names=1)
D3_CHG = read.table("D3.methyl.CHG",row.names=1)
D4_CHG = read.table("D4.methyl.CHG",row.names=1)
D3_CHH = read.table("D3.methyl.CHH",row.names=1)
D4_CHH = read.table("D4.methyl.CHH",row.names=1)


mat = cbind(D3_CG, D4_CG,D3_CHG, D4_CHG,D3_CHH, D4_CHH)
colnames(mat) =c("Meth","D3_CG","D4_CG","D3_CHG","D4_CHG","D3_CHH","D4_CHH")
m <- melt(mat,id=c("Meth"))

colnames(m) = c("Methylation_levels", "Sample", "Distribution")
m$Distribution = log2(m$Distribution)

p <- ggplot(m, aes(x=Methylation_levels,y=Distribution, col=Sample)) +  geom_line(size=1.2, alpha=0.5)
#p <- ggplot(m, aes(x=Methylation_levels,y=Distribution, fill=Sample)) +  geom_bar(stat="identity",alpha=0.5,position="dodge")
p +  theme( text = element_text(size=30),
            axis.line = element_line(size = 1),
            panel.background = element_blank()) +
  labs(title="Distribution of C site Methylation Levels ", x="Methylation Levels", y="log2(Counts)")


## nature 2010
setwd('/lustre/user/liclab/zhangc/proj/lingte/RRBS/STAT')
D3_CG = read.table("D3.methy.CG.percent")
D4_CG = read.table("D4.methy.CG.percent",row.names=1)
D3_CG$V2 = D3_CG$V2/sum(D3_CG$V2) * 100
D4_CG$V2 = D4_CG$V2/sum(D4_CG$V2) * 100
mat = cbind(D3_CG, D4_CG)
colnames(mat) = c("Methylation","D3", "D4")
mat.m = melt(mat, id=c("Methylation"))

p <- ggplot(mat.m, aes(x=Methylation, y=value, fill=variable)) +  geom_bar(stat="identity",position="dodge")
p +  theme( text = element_text(size=30),
            axis.line = element_line(size = 1),
            panel.background = element_blank()) +
  labs( x="Methylation (%)", y="CpGs (%)")



mat = as.matrix(t(mat[,c("D3","D4")]))
par(cex.axis=.4, mgp=c(2,0.5,0),tck=-0.02, font.axis=2, cex.lab=0.5, font.lab=2)
barplot(mat,col=c("red","darkblue"), beside=TRUE,axis.lty=1 ,width=rep(1.5,42),space=c(0,0.5),border=NA,xlim=c(0,80),ylim=c(0,40),xlab="Methylation (%)",ylab="CpGs (%)")




###################  DMR gene region pie plot

setwd("/lustre/user/liclab/zhangc/proj/lingte/RRBS/DMR/bismark/CG_w100p5_f1d0.2/filter/0.3")

# DMR
df <- data.frame(type=c("upstream","downstream","exonic","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","UTR3","UTR5"),
                 numbers=c(185, 53, 206, 529, 564, 30, 74, 33,50) )

# DMC
df <- data.frame(type=c("upstream","downstream","exonic","intergenic","intronic","ncRNA_exonic","ncRNA_intronic","UTR3","UTR5"),
                 numbers=c(142, 69, 314, 1117, 1030, 48, 103, 45,61) )


library(ggplot2)

p = ggplot(df, aes(x="", y=numbers, fill=type)) +
  geom_bar(width=1,stat="identity") +
  coord_polar("y",start=0) +
  theme( text = element_text(size=20),
         panel.background=element_blank(),
         panel.border = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks = element_blank(),
         panel.grid.major = element_line(colour = "gray")
  )+
  ggtitle("Proportion of DMCs distribution")



############  violin plot
setwd("/lustre/user/liclab/zhangc/proj/lingte/RRBS/DMR/bismark/CG_w100p5_f1d0.2/filter/0.3/plot")
D3 = read.table("D3.violin")$V3
D4 = read.table("D4.violin")$V3

mat = cbind(D3,D4)

mat.melt = melt(mat)

colnames(mat.melt) = c("Num", "Sample", "Methylation_Level")
p = ggplot(mat.melt, aes(Sample, Methylation_Level)) + geom_violin(aes(fill = Sample)) +
  theme( text = element_text(size=25),
         panel.background=element_blank(),
         axis.line = element_line(size = 1, colour="gray")) +
  labs(title="Methylation levels in DMRs\n(P < 2.2e-16) ",x="Samples", y="Methylation Levels")

par(mar=c(5.1,5.1,4.1,2.1),cex=1.5)
library(vioplot)
plot(0,0,xlim=c(0,3),ylim=c(0,1),main="Methylation levels in DMRs",xlab="(P < 2.2e-16)",ylab="Methylation level",col="white", xaxt="n")
vioplot(D3,at=1, col="gold", add=T)
vioplot(D4,at=2, col="gold", add=T)
axis(1,at=1:2, labels=c("D3","D4"))

##########  DMC

df <- data.frame(type=c("CG","CHG","CHH"),
                 numbers=c(2616, 94, 207) )

p = ggplot(df, aes(x="", y=numbers, fill=type)) +
  geom_bar(width=1,stat="identity") +
  coord_polar("y",start=0) +
  theme( text = element_text(size=20),
         panel.background=element_blank(),
         panel.border = element_blank(),
         axis.title.x = element_blank(),
         axis.title.y = element_blank(),
         axis.ticks = element_blank(),
         panel.grid.major = element_line(colour = "gray")
  )+
  ggtitle("Counts of DMC")


