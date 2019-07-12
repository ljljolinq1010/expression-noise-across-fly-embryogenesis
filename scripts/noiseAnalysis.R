library(edgeR) 
library(sva) 
library(data.table) 
library(gridExtra) 
library(dendextend) 
library(lawstat) 
library(RColorBrewer) 
library(gplots) 
library(pheatmap) 
library(ppcor) 
library(plyr) 
library(ggplot2) 
library(gdata) 
## plot color
pal <- c(brewer.pal(9, "Set1"), brewer.pal(8, "Set2"), brewer.pal(12, "Set3"))
#####** data input **#####
## Load count summary
count.E1_1<-read.csv("data/BRBseq/E1_1_count.txt", sep="\t", header=TRUE, row.names=1, check.names = F)
count.E1_2<-read.csv("data/BRBseq/E1_2_count.txt", sep="\t", header=TRUE, row.names=1, check.names = F)
count.E2<-read.csv("data/BRBseq/E2_count.txt", sep="\t", header=TRUE, row.names=1, check.names = F)
count.E3<-read.csv("data/BRBseq/E3_count.txt", sep="\t", header=TRUE, row.names=1, check.names = F)
count.all<-cbind(count.E1_1,count.E1_2,count.E2,count.E3)
## order sample based on development stages
stage.name<-c("embryo2-3h","embryo5-6h","embryo8-9h","embryo11-12h","embryo14-15h","embryo17-18h","embryo20-21h","embryo23-24h")
library.name<-c("E1_1","E1_2","E2","E3")
stage.id<-paste0(rep("E",8),c(1:8))
count.ordered<-count.all[,FALSE]
for (i in c(1:length(stage.name))){
  count.ordered.temp<-count.all[,grepl(stage.name[i],colnames(count.all))]
  count.ordered<-cbind(count.ordered,count.ordered.temp)
}
count.all<-count.ordered
## only use protein coding genes
gene.type<-read.csv("data/annotation_files/GeneID_Type.91.txt", sep="\t", header=TRUE, row.names=1, check.names = F)
count.all<-count.all[which(gene.type=="protein_coding"),]
count.all<-na.omit(count.all)

#####** process expression data **#####
#####* quality check *#####
## sample filter
reads.numb = colSums(count.all)/1E6
expressed.gene.numb = colSums(count.all[,names(count.all)] > 0)
sample.quality<-data.frame(reads.numb,expressed.gene.numb)
sample.name.retained<-names(reads.numb[reads.numb>0.3&expressed.gene.numb>4500])
count.retained<-count.all[,sample.name.retained]

## plot the number of reads and of expressed genes
par(mfrow=c(1,1))
par(mar=c(7,5,2,2))
plot(sample.quality,xlim=c(0,18) ,pch=16,ylim=c(0,max(expressed.gene.numb)), xlab = "Uniquely mapped reads number [millions]", ylab = "Expressed genes number [count >= 1]")
# mark outliers
sample.id.retained=which(sample.quality$reads.numb>0.3 & sample.quality$expressed.gene.numb>4500)
sample.id.filtered<-which(sample.quality$reads.numb<=0.3 | sample.quality$expressed.gene.numb<=4500)
points( sample.quality[sample.id.retained,], col="orange",pch=19 )

## gene filter
# use cpm to remove lowly expressed genes (mean cpm across samples >1)
count.retained.dge = DGEList(counts=count.retained)
count.retained.dge = count.retained.dge[rowMeans(cpm(count.retained.dge))>1, ]
# keep the original count value, not the cpm value
count.retained = count.retained[rownames(count.retained.dge$counts), ]

#####* normalization *#####
## Voom normalization (quantile) 
count.retained.norm = voom(counts = count.retained, normalize.method="quantile", plot=F)
count.retained.norm<-count.retained.norm$E
## set negative value as 0
count.retained.norm[count.retained.norm<0]<-0

#####* remove potential batch effect (combat) *#####
## batch information (four libraries)
batch.numb<-list(NA,NA,NA,NA,NA,NA,NA,NA)

for (i in 1:length(stage.name)) {
  for (j in 1:length(library.name)) {
    batch.numb[[i]][j]<-ncol(data.frame(count.retained.norm[,grepl(paste0(stage.name[i],"_",library.name[j]),colnames(count.retained.norm))]))
  }
}
batch.file<-unlist(lapply(batch.numb,function(x) c(rep("A",x[1]),rep("B",x[2]),rep("C",x[3]),rep("D",x[4]))))
## parametric adjustment
combat.retained<-ComBat(dat=count.retained.norm, batch=batch.file, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
## after parametric adjustment, if some value <0, change it into 0
combat.retained[combat.retained<0]<-0
## for genes with value equal to 0 before parametric adjustment, change it into 0
combat.retained[which(count.retained.norm==0)]<-0

#####** data visualization **#####
data.input.retained=data.frame(combat.retained)
stage.color<-c("red" ,"orange" , "green","cyan","blue","purple","forestgreen","grey")
replicates.numb.retained<-c(34, 34, 36, 25, 26, 22 ,30 ,32)
colnames(data.input.retained)<-rep(stage.id,times=replicates.numb.retained)
cor.dist.matrix.retained = as.dist((1 - cor(data.input.retained))/2)
fit.retained <- cmdscale(cor.dist.matrix.retained, eig=TRUE, k=2) # k is the number of dim
par( mar = c(5,5,4,2))
plot(fit.retained$points[,1:2],col=rep(stage.color,times=replicates.numb.retained), 
     xlab="Dim 1", ylab="Dim 2",pch=17, main="MDS",cex=1.5,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,xlim=c(-0.12,0.25))
legend("topright",legend=stage.id, 
       col=stage.color,pch=rep(17,times=8),cex=1.5,bty="n")

#####** test why two clusters **#####
## minor dataset (the samples without developmental trajectory)
minor.dataset<-data.input.retained[fit.retained$points[,1]<0 & fit.retained$points[,2]>0]
minor.dataset$Ensembl.Gene.ID<-row.names(minor.dataset)
minor.dataset<-minor.dataset[c(90,1:89)]
## major dataset (the samples with good developmental trajectory)
major.dataset<-data.input.retained[!(fit.retained$points[,1]<0 & fit.retained$points[,2]>0)]
major.dataset$Ensembl.Gene.ID<-row.names(major.dataset)
major.dataset<-major.dataset[c(151,1:150)]

#####* expression correlation with unfertilized egg *#####
exp.egg<-read.table("data/unfertilized_egg/GSE68062_Gene_abundances_after_FPKM_normalization.txt",sep = "\t",h=T)
exp.egg<-exp.egg[c("mel","mel.E")]
colnames(exp.egg)<-c("Ensembl.Gene.ID","egg")
## minor dataset with unfertilized egg 
minor.dataset.merged<-merge(minor.dataset,exp.egg,by="Ensembl.Gene.ID")
minor.cor<-list()
for (i in c(1:8)) { 
  minor.cor[[i]]<-apply(minor.dataset.merged[,grepl(stage.id[i],colnames(minor.dataset.merged))],2,function(x) cor(x,minor.dataset.merged$egg,method = "spearman"))
}

## major dataset with unfertilized egg 
major.dataset.merged<-merge(major.dataset,exp.egg,by="Ensembl.Gene.ID")
major.cor<-list()
for (i in c(1:8)) { 
  major.cor[[i]]<-apply(major.dataset.merged[,grepl(stage.id[i],colnames(major.dataset.merged))],2,function(x) cor(x,major.dataset.merged$egg,method = "spearman"))
}
## plot
par(mar=c(7,5,2,2))
boxplot(minor.cor, xaxt = "n", col=pal[2], boxwex=0.35, ylim=c(0.2,0.85),cex.main=1.5,cex.axis=1.5,
        cex.lab=1.5,pch=16,outcex=0.35,ylab="Spearman's Rho",at = 1:8 - 0.18) 
boxplot(major.cor, xaxt = "n", add = TRUE, col=pal[1], boxwex=0.35,pch=16,outcex=0.35,axes=FALSE,
        at = 1:8 + 0.18) 
text(x =1:8+0.18, y = 0.15, srt = 45,cex=1.5, adj = 1,  labels = stage.id,xpd = TRUE)

#####* compare the expression of meiosis related genes *#####
## meiosis genes
meiosisI.genes<-read.table("data/annotation_files/meiosis_I.txt")
meiosisII.genes<-read.table("data/annotation_files/meiosis_II.txt")
meiosis.genes<-rbind(meiosisI.genes,meiosisII.genes)
names(meiosis.genes)<-"Ensembl.Gene.ID"
meiosis.genes.exp.minor<-merge(meiosis.genes,minor.dataset,by="Ensembl.Gene.ID")
colnames(meiosis.genes.exp.minor)[-1]<-c(1:89)
meiosis.genes.exp.major<-merge(meiosis.genes,major.dataset,by="Ensembl.Gene.ID")
colnames(meiosis.genes.exp.major)[-1]<-c(90:239)

## pheatmap, the legend is scaled value, (x-mean)/sd
sample.annotation <- data.frame(c(rep(stage.id,times=c(23,10,9,11,5,4,9,18)),rep(stage.id,times=c(11,24,27,14,21,18,21,14))))
colnames(sample.annotation)<-"Stage"
rownames(sample.annotation) <- c(1:239) 
Stage        <- rep(stage.color,times=1)
names(Stage) <- rep(stage.id,times=1)
anno_colors <- list(Stage = Stage)
pheatmap(as.matrix(cbind(meiosis.genes.exp.minor[-1],meiosis.genes.exp.major[-1])), cluster_cols=F,cluster_rows=F,scale='row',main="Egg activation related genes",show_rownames=F,
         fontsize=20,labels_col = "",color=rev(brewer.pal(11,"RdBu")),annotation = sample.annotation,annotation_colors = anno_colors)

#####** expression noise across development **#####
#####* adjusted SD *#####
##### all genes after quality control #####
noise.ajsd <-noiseAJSD(major.dataset[-1],stage.id)
## plot
boxplot(noise.ajsd,outline=FALSE,ylim=c(-0.1,2),notch=T,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="Global adjusted SD",
        col=stage.color, xaxt = "n")

##### only use genes expressed in all stages #####
noise.ajsd.df<-data.frame(rownames(major.dataset[-1]),noise.ajsd)
rownames(noise.ajsd.df)<-NULL
names(noise.ajsd.df)<-c("Ensembl.Gene.ID",stage.id)

## mean value in each stage
mean.exp.list<-list()
for (i in c(1:8)) {
  sub.data<-major.dataset[,grepl(stage.id[i],colnames(major.dataset))]
  mean.exp.list[[i]]<-rowMeans(sub.data)
}
major.dataset.mean<-data.frame(rownames(major.dataset),mean.exp.list)
names(major.dataset.mean)<-c("Ensembl.Gene.ID",stage.id)
rownames(major.dataset.mean)<-NULL

## define expressed genes (expression value > 1)
major.dataset.mean$freq<-apply(major.dataset.mean[,c(2:9)],1,function(x) length(x[x>1]))
expGene<-subset(major.dataset.mean,major.dataset.mean$freq>7)
expGeneNoise<-noise.ajsd.df[noise.ajsd.df$Ensembl.Gene.ID%in%expGene$Ensembl.Gene.ID,]

## plot
boxplot(expGeneNoise[-1],outline=FALSE,ylim=c(-0.1,2),notch=T,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="Global adjusted SD",
        col=stage.color, xaxt = "n")

##### only use genes with constant expression across development #####
sigValue<-c()
for (i in c(1:8004)) {
  a<-major.dataset[i,-1]
  b<-data.frame(as.numeric(a),c(rep(stage.id,times=c(11,24,27,14,21,18,21,14))))
  colnames(b)<-c("exp","stage")
  # compute the analysis of variance
  res.aov <- aov(exp ~ stage, data = b)
  # summary of the analysis
  results<-summary(res.aov)
  sigValue[i]<-results[[1]][1,5]
}
sigValue<-data.frame(sigValue)
## qvalue
sigValueCorrect<- qvalue(sigValue$sigValue,pi0=1)$qvalues
sigValueCorrect<-data.frame(rownames(major.dataset),as.numeric(sigValueCorrect))

## constant expression genes defined as qvalue >0.05
constantGene<-subset(sigValueCorrect,sigValueCorrect$as.numeric.sigValueCorrect.>0.05)
constantGeneNoise<-noise.ajsd.df[noise.ajsd.df$Ensembl.Gene.ID%in%constantGene$rownames.major.dataset.,]
## plot
boxplot(constantGeneNoise[-1],outline=FALSE,ylim=c(-0.1,2),notch=T,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="Global adjusted SD",
        col=stage.color, xaxt = "n")

##### transcription factor #####
TF<-read.table("data/annotation_files/tf.txt")
TFnoise<-noise.ajsd.df[noise.ajsd.df$Ensembl.Gene.ID%in%TF$V1,]
## plot
boxplot(TFnoise[-1],outline=FALSE,ylim=c(-0.1,2),notch=T,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="Global adjusted SD",
        col=stage.color, xaxt = "n")

##### broad promter and core promoter noise  #####
promoterSI<-fread("data/annotation_files/promoter_si.txt")
colnames(promoterSI)<-c("Gene.Name","SI")
gene.id.name<-read.table("data/annotation_files/fly_ensemblID_geneName.txt",sep="\t",quote="",h=T)
names(gene.id.name)<-c("Ensembl.Gene.ID","Gene.Name")
promoterSI<-merge(promoterSI,gene.id.name,by="Gene.Name")
promoterSInoise<-merge(noise.ajsd.df,promoterSI,by="Ensembl.Gene.ID")

broadPromoter<-subset(promoterSInoise,promoterSInoise$SI<=-1.5)
corePromoter<-subset(promoterSInoise,promoterSInoise$SI>-1.5)

## plot both together
noiseData<-list(broadPromoter[,2],corePromoter[,2],broadPromoter[,3],corePromoter[,3],broadPromoter[,4],corePromoter[,4],
                broadPromoter[,5],corePromoter[,5],broadPromoter[,6],corePromoter[,6],broadPromoter[,7],corePromoter[,7],
                broadPromoter[,8],corePromoter[,8],broadPromoter[,9],corePromoter[,9])
noisePlot<-boxplot(noiseData,at=c(1,2,4,5,7,8,10,11,13,14,16,17,19,20,22,23),outline=FALSE,ylim=c(-0.1,2),notch=T,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="Adjusted SD",
                   col=rep(stage.color,each=2), xaxt = "n")


#####* distance to median *#####
noise.dm <- noiseDM(major.dataset[-1],stage.id)
## plot
boxplot(noise.dm,outline=FALSE,ylim=c(-2.2,2),notch=T,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="Distance to Median (DM)",
        col=stage.color, xaxt = "n")

#####* cv *#####
noise.cv<-list()
for (i in 1:length(stage.id)) {
  subData<-major.dataset[-1][,grepl(stage.id[i],colnames(major.dataset[-1]))]
  subData$mean<-rowMeans(subData)
  noise.cv[[i]] <-apply(subData[,-(ncol(subData))],1,function(x) sd(x)/mean(x))
}
noise.cv.df<-data.frame(rownames(major.dataset),noise.cv)
rownames(noise.cv.df)<-NULL
names(noise.cv.df)<-c("Ensembl.Gene.ID",stage.id)
## plot
boxplot(noise.cv,outline=FALSE,ylim=c(0,3),notch=T,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="CV",
        col=stage.color, xaxt = "n")


#####** why lower noise in phylotypic stage E3? (histone modifications) **#####
#####* histone modification signal (Z score) across development *#####
promoter.hist<-fread("data/histone_modifications/promoter_histone_zscore.wig",sep="\t")
colnames(promoter.hist)<-c("chr","start","end","Ensembl.Gene.ID",
                           paste0(rep( c("pro_H3K27Ac_","pro_H3K4Me1_","pro_H3K4Me3_",
                                         "pro_H3K9Ac_"),each=6), rep(c(1,2,3,4,5,6),4)))
promoter.hist<-data.frame(promoter.hist)
promoter.hist<-promoter.hist[promoter.hist$Ensembl.Gene.ID%in%noise.ajsd.df$Ensembl.Gene.ID,]
for  (i in 5:28) {
  promoter.hist[,i]<-as.numeric(promoter.hist[,i])
}
for  (i in 5:28) {
  promoter.hist<-promoter.hist[!is.na(promoter.hist[,i]),]
}

## gene body signal
gene.hist<-fread("data/histone_modifications/genebody_histone_zscore.wig",sep="\t")
colnames(gene.hist)<-c("chr","start","end","Ensembl.Gene.ID",
                       paste0(rep( c("gene_H3K27Ac_","gene_H3K4Me1_","gene_H3K4Me3_",
                                     "gene_H3K9Ac_"),each=6), rep(c(1,2,3,4,5,6),4)))
gene.hist<-data.frame(gene.hist)
for  (i in 5:28) {
  gene.hist[,i]<-as.numeric(gene.hist[,i])
}
for  (i in 5:28) {
  gene.hist<-gene.hist[!is.na(gene.hist[,i]),]
}
gene.hist<-gene.hist[gene.hist$Ensembl.Gene.ID%in%noise.ajsd.df$Ensembl.Gene.ID,]
## plot: gene body signal for H3K9Ac
boxplot(gene.hist[,c(23:28)],notch=T,outline=FALSE,pch=16,outcex=0.5,boxwex=0.7, 
        col=pal[2],xaxt = "n",ylim=c(-5,15),ylab="Zscore (gene body)",main="H3K9Ac",cex.axis = 2,cex.main = 2,cex.lab = 2)

#####* histone modification signal (tag)  and noise correlation *#####
## promoter signal 
promoter.hist<-fread("data/histone_modifications/promoter_histone_tag.wig",sep="\t")
colnames(promoter.hist)<-c("chr","start","end","Ensembl.Gene.ID",
                           paste0(rep( c("pro_H3K27Ac_","pro_H3K4Me1_","pro_H3K4Me3_",
                                         "pro_H3K9Ac_"),each=6), rep(c(1,2,3,4,5,6),4)))
promoter.hist<-data.frame(promoter.hist)
for  (i in 5:28) {
  promoter.hist[,i]<-as.numeric(promoter.hist[,i])
}
for  (i in 5:28) {
  promoter.hist<-promoter.hist[!is.na(promoter.hist[,i]),]
}

## gene body signal
gene.hist<-fread("data/histone_modifications/genebody_histone_tag.wig",sep="\t")
colnames(gene.hist)<-c("chr","start","end","Ensembl.Gene.ID",
                       paste0(rep( c("gene_H3K27Ac_","gene_H3K4Me1_","gene_H3K4Me3_",
                                     "gene_H3K9Ac_"),each=6), rep(c(1,2,3,4,5,6),4)))
gene.hist<-data.frame(gene.hist)
for  (i in 5:28) {
  gene.hist[,i]<-as.numeric(gene.hist[,i])
}
for  (i in 5:28) {
  gene.hist<-gene.hist[!is.na(gene.hist[,i]),]
}

## merge histone and noise  
promoter.hist.noise<-merge(promoter.hist,noise.ajsd.df,by="Ensembl.Gene.ID")
gene.hist.noise<-merge(gene.hist,noise.ajsd.df,by="Ensembl.Gene.ID")

## mean signal and mean noise correlation
histone.markers<-c("H3K4Me1","H3K4Me3","H3K9Ac","H3K27Ac")
# promoter
promoter.hist.mean<-list()
for (i in c(1:4)) {
  promoter.temp.data<-promoter.hist.noise[,paste0("pro_",rep(histone.markers[i],6),"_",c(1:6))]
  promoter.hist.mean[[i]]<-rowMeans(promoter.temp.data)
}
promoter.noise.mean<-rowMeans(promoter.hist.noise[,stage.id])
promoter.hist.noise.mean<-data.frame(promoter.hist.noise$Ensembl.Gene.ID,data.frame(promoter.hist.mean),data.frame(promoter.noise.mean))
colnames(promoter.hist.noise.mean)<-c("Ensembl.Gene.ID",histone.markers,"noise")

# gene body
gene.hist.mean<-list()
for (i in c(1:4)) {
  gene.temp.data<-gene.hist.noise[,paste0("gene_",rep(histone.markers[i],6),"_",c(1:6))]
  gene.hist.mean[[i]]<-rowMeans(gene.temp.data)
}
gene.noise.mean<-rowMeans(gene.hist.noise[,stage.id])
gene.hist.noise.mean<-data.frame(gene.hist.noise$Ensembl.Gene.ID, data.frame(gene.hist.mean),data.frame(gene.noise.mean))
colnames(gene.hist.noise.mean)<-c("Ensembl.Gene.ID",histone.markers,"noise")

## correlation 
cor.results<-mean_cor(promoter.hist.noise.mean[,c(2:6)],gene.hist.noise.mean[,c(2:6)],histone.markers)

## plot
promoter.cor.eff<-cor.results[[1]]
gene.cor.eff<-cor.results[[3]]
signif.marks<-c(rep("***",8))
cor.eff<-c(promoter.cor.eff[1],gene.cor.eff[1],promoter.cor.eff[2],gene.cor.eff[2],
           promoter.cor.eff[3],gene.cor.eff[3],promoter.cor.eff[4],gene.cor.eff[4])

cor.eff.plot<-barplot(cor.eff,main="Correlation with promoter Phastcons", horiz=TRUE,xlim=c(-1,0.1),
                      names.arg =c("H3K4Me1","","H3K4Me3","","H3K9Ac","","H3K27Ac",""),xlab = "Spearman's Rho",
                      cex.axis =1.3,cex.names = 1.3,cex.lab=1.3,width=c(3,3,3,3),col=rep(c(pal[1],pal[2]),4))

legend("topleft",c("promoter","gene body"), fill=c(pal[1],pal[2]), horiz=TRUE, cex=1,bty = "n")
for (i in c(1:length(signif.marks))) {
  text(0.05,cor.eff.plot[i],signif.marks[i],cex = 1.2) ## * <0.05, ** <0.01, *** <0.001
}

#####** histone and promoter sequence conservation correlation **#####
## promoter phastCons
gene.id.name<-read.table("data/promoter/fly_ensemblID_geneName.txt",sep="\t",quote="",h=T)
names(gene.id.name)<-c("Ensembl.Gene.ID","Gene.Name")
promoter.phastCons<-fread("data/promoter/promoter_narrowPhast.txt")

names(promoter.phastCons)<-c("Gene.Name","PhastCons")
promoter.phastCons<-promoter.phastCons[!grepl("^.$", promoter.phastCons$PhastCons),]
promoter.phastCons$PhastCons<-as.numeric(promoter.phastCons$PhastCons)
promoter.phastCons<-merge(promoter.phastCons,gene.id.name,by="Gene.Name")
promoter.phastCons$Gene.Name<-NULL
## merge
promoter.hist.noise.mean.phastCons<-merge(promoter.hist.noise.mean, promoter.phastCons,by="Ensembl.Gene.ID")
promoter.hist.noise.mean.phastCons<-na.omit(promoter.hist.noise.mean.phastCons)
gene.hist.noise.mean.phastCons<-merge(gene.hist.noise.mean, promoter.phastCons,by="Ensembl.Gene.ID")
gene.hist.noise.mean.phastCons<-na.omit(gene.hist.noise.mean.phastCons)

## correlation
cor.results<-mean_cor(promoter.hist.noise.mean.phastCons[,c(2:5,7)],gene.hist.noise.mean.phastCons[,c(2:5,7)],histone.markers)


#####* histone and promoter sequence pi correlation *#####
## the pi caluclated by using all snps or only common snps (remove minor allele frequence <0.05) is quite similar
promoter.pi<-fread("data/promoter/promoter_narrow_pi.txt")
promoter.pi<-promoter.pi[,c(4,5)]
names(promoter.pi)<-c("Gene.Name","pi")
promoter.pi<-promoter.pi[!grepl("^.$", promoter.pi$pi),]
promoter.pi$pi<-as.numeric(promoter.pi$pi)
promoter.pi<-merge(promoter.pi,gene.id.name,by="Gene.Name")
promoter.pi$Gene.Name<-NULL
promoter.pi<-na.omit(promoter.pi)
promoter.pi$pi<-promoter.pi$pi*100000/60

## merge
promoter.hist.noise.mean.pi<-merge(promoter.hist.noise.mean, promoter.pi,by="Ensembl.Gene.ID")
promoter.hist.noise.mean.pi<-na.omit(promoter.hist.noise.mean.pi)
gene.hist.noise.mean.pi<-merge(gene.hist.noise.mean, promoter.pi,by="Ensembl.Gene.ID")
gene.hist.noise.mean.pi<-na.omit(gene.hist.noise.mean.pi)
## correlation
cor.results<-mean_cor(promoter.hist.noise.mean.pi[,c(2:5,7)],gene.hist.noise.mean.pi[,c(2:5,7)],histone.markers)

#####** compare sequence conservation of promoters of stage specifically expressed genes **#####
load("data/expression_modules/stage_specific_high_exp_genes.RData")

## get phastcons for stage specifically expressed genes
cluster.phastCons<-list()
for (i in 1:8) {
  temp.data<-exp.clusters[[i]][1]
  names(temp.data)<-"Ensembl.Gene.ID"
  temp.data<-promoter.phastCons$PhastCons[promoter.phastCons$Ensembl.Gene.ID%in%temp.data$Ensembl.Gene.ID]
  names(temp.data)<-"PhastCons"
  cluster.phastCons[[i]]<-temp.data
}

## plot
par(mar=c(7,5,2,2))
phast.plot<-boxplot(cluster.phastCons,outline=FALSE,ylim=c(-0.1,1.2),notch=T,cex.main=1.5,cex.axis=1.5,cex.lab=1.5,ylab="PhasctCons score",
                    col=stage.color,xaxt = "n")
text(x =phast.plot$names, y = -0.2, srt = 45,cex=1.5, adj = 1,  labels = stage.id,xpd = TRUE)

