#####* function used to calculate noise: global adjusted standard deviation *#####
predictSD<-function(expData) {
  ## polynomial model 
  expData$mean<-rowMeans(expData)
  expData$sd<-apply(expData[,-ncol(expData)],1,function(x) sd(x))
  m1 <- lm(sd~mean, expData)
  m2 <- update(m1, .~. + I(mean^2), expData)
  m3 <- update(m2, .~. + I(mean^3), expData)
  m4 <- update(m3, .~. + I(mean^4), expData)
  m5 <- update(m4, .~. + I(mean^5), expData)
  m6 <- update(m4, .~. + I(mean^6), expData)
  totalM<-list(m1,m2,m3,m4,m5,m6)
  for (i in c(1:5)) {
    modelComp <- anova(totalM[[i]],totalM[[i+1]])
    if (is.na(modelComp$`Pr(>F)`[2]) | modelComp$`Pr(>F)`[2]>0.05) {
      m<-totalM[[i]]
      break
    }
  }
  return(predict(m))
}

noiseAJSD<-function(expData,stageID) {
  expectSD<-predictSD(expData)
  ajsdList<-list()
  for (i in 1:length(stageID)) {
    subData<-expData[,grepl(stageID[i],colnames(expData))]
    subData$mean<-rowMeans(subData)
    subData$sd<-apply(subData[,-(ncol(subData))],1,function(x) sd(x))
    ajsdList[[i]]<-subData$sd/expectSD
  }  
  return(ajsdList)
}

#####* function used to calculate noise: global adjusted distance to median *#####
library(scran) 
DM <- function(globleMean, globleCVs, subDataCVs, win.size=51) {
  ## compute the distance to median of the CV2 (squared coefficient of variation)
  ## globleMean: mean expression of all samples
  ## globleCVs: CVs of all samples
  ## subDataCVs: CVs of a specific stage
  keep <- globleMean > 0 & !is.na(globleCVs) & globleCVs > 0 & !is.na(subDataCVs) & subDataCVs>0
  globleMean.expr <- globleMean[keep]
  globleCVs.expr <- log10(globleCVs[keep])
  subDataCVs.expr <- log10(subDataCVs[keep])
  o <- order(globleMean.expr)
  if (win.size%%2L==0L) {
    win.size <- win.size+1L
  }
  med.trend <- runmed(globleCVs.expr[o], k=win.size)
  med.trend[o] <- med.trend
  dm.out <- subDataCVs.expr - med.trend
  DM <- rep(NA_real_, length(keep))
  DM[keep] <- dm.out
  return(DM)
}

noiseDM<-function(expData,stageID) {
  dmList<-list()
  sampleNumbList<-list()
  globleMean<-rowMeans(expData)
  globleCVs<-apply(expData,1,function(x) (sd(x)/mean(x))^2)
  
  for (i in 1:length(stageID)) {
    subData<-expData[,grepl(stageID[i],colnames(expData))]
    subData$mean<-rowMeans(subData)
    subData$cvs<-apply(subData[,-(ncol(subData))],1,function(x) (sd(x)/mean(x))^2)
    dmList[[i]]<-DM(globleMean,globleCVs,subData$cvs,win.size=51)
  }  
  return(dmList)
}

#####* promoter and gene body histone signal based correlation *#####
mean_cor<-function(promoter.data,gene.data,histone.markers) {
  promoter.cor.eff<-c()
  promoter.cor.pvalue<-c()
  gene.cor.eff<-c()
  gene.cor.pvalue<-c()
  n=length(histone.markers)
  colnames(promoter.data)<-c(histone.markers,"variable")
  colnames(gene.data)<-c(histone.markers,"variable")
  
  for (i in c(1:n)) { ## four markers
    ## promoter
    promoter.cor.results<-cor.test(promoter.data[,histone.markers[i]],promoter.data$variable,method="spearman")
    promoter.cor.eff[i]<-promoter.cor.results$estimate
    promoter.cor.pvalue[i]<-promoter.cor.results$p.value
    ## gene body
    gene.cor.results<-cor.test(gene.data[,histone.markers[i]],gene.data$variable,method="spearman")
    gene.cor.eff[i]<-gene.cor.results$estimate
    gene.cor.pvalue[i]<-gene.cor.results$p.value
  }
  results<-list(promoter.eff=promoter.cor.eff,promoter.pvalue=promoter.cor.pvalue,
                gene.eff=gene.cor.eff,gene.pvalue=gene.cor.pvalue)
  return(results)
}

