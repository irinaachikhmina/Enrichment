#Integration of omics data
#Irina Makarova


library(iClusterPlus)
library(ggfortify)
library(parallel)

#Скачиваем данные
mrna_d = read.table('mrna.tsv', sep = '\t', header = T)
mrna = data.frame(mrna_d[-1])
row.names(mrna) = mrna_d[['Gene']]
mrna = t(mrna)

mirna_d = read.table('mirna.tsv', sep = '\t', header = T)
mirna = data.frame(mirna_d[-1])
row.names(mirna) = mirna_d[['miRNA']]
mirna = t(mirna)

prot_d = read.table('prot.tsv', sep = '\t', header = T)
prot = data.frame(prot_d[-1])
row.names(prot) = prot_d[['Protein']]
prot = t(prot)

#Предварительный анализ данных
##mRNA
#####PCA по образцам
modelPCA = prcomp(mrna)
plot(modelPCA)
summary(modelPCA)
modelPCA$rotation
autoplot(modelPCA, data = mrna, label = TRUE, label.size = 4)
#####k-means по образцам
autoplot(kmeans(mrna, 3), data = mrna, label = TRUE, label.size = 3) #пробуем 3 кластера, один из них неплохо выделяется, есть образец-потенциальный выброс (А133, аналогично PCA).
autoplot(kmeans(mrna, 2), data = mrna, label = TRUE, label.size = 3) #образцы лучше кластеризиуются на 2 группы, чем на 3
autoplot(kmeans(mrna, 4), data = mrna, label = TRUE, label.size = 3)
#####boxplot в разрезе образцов
boxplot(t(mrna)) #A133
####Гены
boxplot(mrna)
#####корреляция
correlation = cor(mrna,method = 'spearman')
heatmap(correlation,
        symm = TRUE,
        distfun = function(x){as.dist(1-x)}) #есть сильно скоррелированные гены, возможно, стоит их кластеризовать и оставить по 1 из кластера

w = which(abs(correlation)>0.8 & row(correlation)<col(correlation), arr.ind=TRUE)
corr_mrna = matrix(colnames(correlation)[w],ncol=2)
corr_mrna
unique(corr_mrna[,1])
#LRRC25 - Low cancer specificity (https://www.proteinatlas.org/ENSG00000175489-LRRC25/pathology), (Leucine Rich Repeat Containing 25) is a Protein Coding gene. Diseases associated with LRRC25 include Campylobacteriosis. Among its related pathways are Vitamin D receptor pathway.
#to remove from the sample
#NCF4 - Low cancer specificity. NCF4 (Neutrophil Cytosolic Factor 4) is a Protein Coding gene. Diseases associated with NCF4 include Granulomatous Disease, Chronic, Autosomal Recessive, 3 and Chronic Granulomatous Disease. Among its related pathways are Signaling by Rho GTPases and Cellular responses to stimuli. Gene Ontology (GO) annotations related to this gene include protein dimerization activity and phosphatidylinositol-3-phosphate binding. 
#to remove from the sample
mrna_sp = subset(mrna, select = -c(ZEB1, HTRA1, FLI1, TCF4, ASPM, CSF1R, APBB1IP, LRRC25, NCF4, COL6A1,  POSTN))
mrna_sp

###Nan нет

##miRNA
####Образцы
modelPCA = prcomp(mirna)
autoplot(modelPCA, data = mirna, label = TRUE, label.size = 3)
autoplot(kmeans(mirna, 3), data = mirna, label = TRUE, label.size = 3) #лучше всего разделяются на 3 кластера
autoplot(kmeans(mirna, 2), data = mirna, label = TRUE, label.size = 3) 
autoplot(kmeans(mirna, 4), data = mirna, label = TRUE, label.size = 3) 

boxplot(t(mirna))
####микроРНК
boxplot(mirna) #сильный разброс
#####корреляция
colnames(mirna) = gsub("-", "_", colnames(mirna))
correlation = cor(mirna,method = 'spearman')
heatmap(correlation,
        symm = TRUE,
        distfun = function(x){as.dist(1-x)})
w = which(abs(correlation)>0.8 & row(correlation)<col(correlation), arr.ind=TRUE)
corr_mirna = matrix(colnames(correlation)[w],ncol=2)
corr_mirna = corr_mirna[-c(27, 38),] #hsa-mir-429, hsa-mir-93 коррелируют с микроРНК, которые были удалены на предыдущем шаге
paste(unique(corr_mirna[,2]), collapse = ", ")
mirna_sp = subset(mirna, 
                  select = -c(hsa_let_7a_2, hsa_let_7a_3, hsa_mir_128_2, hsa_mir_134, hsa_mir_155, hsa_mir_16_2, hsa_mir_181b_1, hsa_mir_181d, hsa_mir_194_2, hsa_mir_199a_2, hsa_mir_199b, hsa_mir_19b_2, hsa_mir_200b, hsa_mir_200c, hsa_mir_20a, hsa_mir_222, hsa_mir_27b, hsa_mir_29b_2, hsa_mir_365_2, hsa_mir_379, hsa_mir_381, hsa_mir_409, hsa_mir_451, hsa_mir_501, hsa_mir_532, hsa_mir_660, hsa_mir_9_2, hsa_mir_92a_1, hsa_mir_92a_2, hsa_mir_99a))


###Nan нет

##Proteins
####Образцы
modelPCA = prcomp(prot)
autoplot(modelPCA, data = prot, label = TRUE, label.size = 3) #ok
boxplot(t(prot))
####Белки
boxplot(prot) #неплохо
#####корреляция
correlation = cor(prot,method = 'spearman')
heatmap(correlation,
        symm = TRUE,
        distfun = function(x){as.dist(1-x)}) #в целом, данные не очень явно кластеризованы
w = which(abs(correlation)>0.8 & row(correlation)<col(correlation), arr.ind=TRUE)
corr_prot = matrix(colnames(correlation)[w],ncol=2)
paste(unique(corr_prot[,2]), collapse = ", ")
prot_sp = subset(prot, 
                  select = -c(ACC_pS79, GATA3, S6_pS240_S244, Chk1_pS345, Mre11))

###Nan
sum(is.nan(prot))

#Scaling, normalization
cat(min(mrna), median(mrna), max(mrna)) #[0, 13]
boxplot(t(mrna))
cat(min(mirna), median(mirna), max(mirna)) #[0, 19]
boxplot(t(mirna))
cat(min(prot), median(prot), max(prot)) #медиана на 0, интервал [-6,6]
boxplot(t(prot))


#Вариант 1(-): Масштабировать по строкам (образцы). В этом случае удастся сохранить паттерн экспрессии образцов  
#Вариант 2(+): Масштабировать по столбцам (гены, микроРНК, белки). В этом случае удастся сохранить соотношения между образцами (сходство, различие), что важнее для задачи кластеризации образцов
mrna_ssp = apply(mrna_sp, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
mrna_s = apply(mrna, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
mirna_ssp = apply(mirna_sp, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
mirna_s = apply(mirna, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
prot_ssp = apply(prot_sp, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
prot_s = apply(prot, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

cat(min(mrna_s), median(mrna_s), max(mrna_s))
boxplot(t(mrna_s))
cat(min(mirna_s), median(mirna_s), max(mirna_s))
boxplot(t(mirna_s))
cat(min(prot_s), median(prot_s), max(prot_s))
boxplot(t(prot_s))

#Outliers mRNA
##образцы
###PCA
modelPCA = prcomp(mrna_s)
plot(modelPCA)
autoplot(modelPCA, data = mrna_s, label = TRUE, label.size = 3) #чуть лучше разделяются, образцы A0D0, A133 выглядят как outlier
summary(modelPCA)
biplot(modelPCA)
hist(modelPCA$rotation[, 1])
modelPCA$rotation[, 1][modelPCA$rotation[, 1] < -0.1]
modelPCA$rotation[, 1][modelPCA$rotation[, 1] > 0.1]
hist(modelPCA$rotation[, 2])

###k-means
dev.off()
autoplot(kmeans(mrna_ssp, 4), data = mrna_ssp, label = TRUE, label.size = 3)
autoplot(kmeans(mrna_ssp, 3), data = mrna_ssp, label = TRUE, label.size = 3)
autoplot(kmeans(mrna_ssp, 2), data = mrna_ssp, label = TRUE, label.size = 3) #лучше всего выделяется в 2 кластера
###boxplot
boxplot(t(mrna_ssp)) #нелохо
###гены
boxplot(mrna_ssp)
medians = apply(mrna_ssp, 2, median)
hist(medians, breaks=20)
medians[medians > .7]

#Outliers - miRNA
###PCA
modelPCA = prcomp(mirna_s)
plot(modelPCA)
autoplot(modelPCA, data = mirna_s, label = TRUE, label.size = 3) #образцы A0D1, A0CD выглядят как outlier
###k-means
autoplot(kmeans(mirna_ssp, 4), data = mirna_ssp, label = TRUE, label.size = 3) 
autoplot(kmeans(mirna_ssp, 3), data = mirna_ssp, label = TRUE, label.size = 3) #лучше всего выделяется в 3 кластера
autoplot(kmeans(mirna_ssp, 2), data = mirna_ssp, label = TRUE, label.size = 3)
###boxplot
boxplot(t(mirna_s))
###микроРНК
medians = apply(mirna_s, 2, median)
hist(medians, breaks=20)
medians[medians > 0.72]

#Outliers - proteins
###PCA
modelPCA = prcomp(prot_s)
plot(modelPCA)
autoplot(modelPCA, data = prot_s, label = TRUE, label.size = 3)
###k-means
autoplot(kmeans(prot_ssp, 4), data = prot_ssp, label = TRUE, label.size = 3)
autoplot(kmeans(prot_ssp, 3), data = prot_ssp, label = TRUE, label.size = 3) #неплохо выделяются 3 кластера
autoplot(kmeans(prot_ssp, 2), data = prot_ssp, label = TRUE, label.size = 3)
###boxplot
boxplot(t(prot_s))
###гены
medians = apply(prot_s, 2, median)
hist(medians, breaks=20)


#Кластеризация по 3 омиксам, iClusterPlus, шкалированные данные, все признаки, n.lambda -значение по умолчанию  (попытка 1)
#код ниже написан на основании руководства к iClusterPlus (iManual())
##Подбор k (количество кастеров = K+1)
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=1, #'mc.cores' > 1 не поддерживается под Windows
                             dt1=mrna_s,
                             dt2=mirna_s,
                             dt3=prot_s,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,
                             n.lambda=NULL, #значение по умолчанию
                             scale.lambda=c(1,1,1),
                             maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

#записываем в otput
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
  }
#nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#находим индексы минимальных значений BIC для каждого значения к
minBICid = apply(BIC,2,which.min)
minBICid

#строим кривую % объясненной разницы
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
  }
devRatMinBIC

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
#кривая немного снижается после 3 кластеров. То есть ожидаемое количество кластеров - 3 (к=2). 

k=2
min(BIC[,k]) #-203296.8

best.fit=output[[k]]$fit[[which.min(BIC[,k])]]

#посмотрим основные признаки, которые были использованы для кластеризации
features = alist()
features[[1]] = colnames(mrna_s)
features[[2]] = colnames(mirna_s)
features[[3]] = colnames(prot_s)
sigfeatures=alist()
for(i in 1:3){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
  }
names(sigfeatures)=c("mrna","mirna","prot")
head(sigfeatures[[1]])
head(sigfeatures[[2]]) #character(0)
head(sigfeatures[[3]])
#видно, что данные микроРНК не были использованы для кластеризации. Пробуем меньшее значение lambda для данных микроРНК
#Рекомендации из мануала: indicates that the selected lambda value is too large and it is not in the same scale as those for the copy number and
#gene expression data. To solve this problem, we need to set the scale.lambda (an argument of tune.iClusterPlus) to a value between 0 and 1.

#lambda (1, 0.05, 1)
setwd('C:/Users/Irina Makarova/Documents/ADBM/System biology - fri/hw3/2lambda_mirna')
for(k in 1:4){ #исключаем вариант 6 кластеров для ускорения расчетов, так как уже на первом этапе видно, что оптимальное количество кластеров 2-4
  cv2.fit = tune.iClusterPlus(cpus=1, #'mc.cores' > 1 не поддерживается под Windows
                              dt1=mrna_s,
                              dt2=mirna_s,
                              dt3=prot_s,
                              type=c("gaussian","gaussian","gaussian"),
                              K=k,
                              n.lambda=101, #меньше, чем отимальное значение (185) для ускорения расчетов
                              scale.lambda=c(1,0.05,1),
                              maxiter=20)
  save(cv2.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv2.fit
}
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)
minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
#
k=2
best.fit=output[[k]]$fit[[which.min(BIC[,k])]]
#основные признаки, которые были использованы для кластеризации
features = alist()
features[[1]] = colnames(mrna_s)
features[[2]] = colnames(mirna_s)
features[[3]] = colnames(prot_s)
sigfeatures=alist()
for(i in 1:3){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("mrna","mirna","prot")
head(sigfeatures[[1]]) 
head(sigfeatures[[2]])
head(sigfeatures[[3]]) #character(0)
#теперь не учтены данные белков. пробуем взять меньшее значение лямбды для набора prot.


#lambda (1, 0.05, 0.1)
setwd('C:/Users/Irina Makarova/Documents/ADBM/System biology - fri/hw3/2lambda_1_005_01_nlNULL')
for(k in 1:4){ #исключаем вариант 6 кластеров для ускорения расчетов, так как уже на первом этапе видно, что оптимальное количество кластеров 2-4
  cv2.fit = tune.iClusterPlus(cpus=1, #'mc.cores' > 1 не поддерживается под Windows
                              dt1=mrna_ssp,
                              dt2=mirna_ssp,
                              dt3=prot_ssp,
                              type=c("gaussian","gaussian","gaussian"),
                              K=k,
                              n.lambda=NULL, 
                              scale.lambda=c(1,0.05,0.1),
                              maxiter=20)
  save(cv2.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv2.fit
}
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)
minBICid = apply(BIC,2,which.min)
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
BIC[minBICid[3]] #-185508.5
plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
#
k=3
best.fit=output[[k]]$fit[[which.min(BIC[,k])]]

#основные признаки, которые были использованы для кластеризации
features = alist()
features[[1]] = colnames(mrna_s)
features[[2]] = colnames(mirna_s)
features[[3]] = colnames(prot_s)
sigfeatures=alist()
for(i in 1:3){
  rowsum=apply(abs(best.fit$beta[[i]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  sigfeatures[[i]]=(features[[i]])[which(rowsum>upper)]
}
names(sigfeatures)=c("mrna","mirna","prot")
head(sigfeatures[[1]]) 
head(sigfeatures[[2]])
head(sigfeatures[[3]])

clusters=getClusters(output)
rownames(clusters)=rownames(mrna_s)
colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
clusters

write.table(paste('cl',best.fit$clusters,sep=""), file="clusterTab.txt",sep='\t',quote=F, col.names = F, row.names = F)

# Смотрим можно ли улучшить модель за счет подбора параметра альфа (в качестве метрики используем BIC) 
# Минимальное значение, найденное на предыдущем этапе,- -185508.5
setwd('C:/Users/Irina Makarova/Documents/ADBM/System biology - fri/hw3/alpha')
mod=iClusterPlus(dt1=mrna_s,
                 dt2=mirna_s,
                 dt3=prot_s,
                 type=c("gaussian","gaussian","gaussian"),
                 #alpha=c(0.5,0.5,0.5),
                 #alpha=c(0.8,0.5,0.5),
                 #alpha=c(0.8,0.5,0.6),
                 alpha=c(0.8,0.1,0.6),
                 lambda=c(1,0.05,0.1),
                 K=3,
                 maxiter=20)

mod$BIC #-202015.9

clusters = paste('cl',mod$clusters,sep="")
write.table(clusters, file="clusterTab.txt",sep='\t',quote=F,row.names=F,col.names = F)
write.table(mod$clusters, file="clusterMembership.txt",sep='\t',quote=F,row.names=F,col.names = F)

gc()

#Кластеризация по 1 омиксу, mrna, iClusterPlus
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=1, 
                             dt1=mrna_s,
                             dt2=NULL,
                             dt3=NULL,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,
                             n.lambda=NULL, #number of points to sample
                             scale.lambda=c(1,1,1),
                             maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

#записываем в otput
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#находим индексы минимальных значений BIC для каждого значения к
minBICid = apply(BIC,2,which.min)
minBICid
k = 2
min(BIC[,k])

#строим кривую % объясненной разницы
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
devRatMinBIC

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")

#Кластеризация по miRNA, iClusterPlus
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=1, 
                             dt1=mirna_s,
                             dt2=NULL,
                             dt3=NULL,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,
                             n.lambda=NULL, #number of points to sample
                             scale.lambda=c(1,1,1),
                             maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

#записываем в otput
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#находим индексы минимальных значений BIC для каждого значения к
minBICid = apply(BIC,2,which.min)
minBICid
k = 2
min(BIC[,k])

#строим кривую % объясненной разницы
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
devRatMinBIC

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")

#Кластеризация по prot, iClusterPlus
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=1, 
                             dt1=prot_s,
                             dt2=NULL,
                             dt3=NULL,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,
                             n.lambda=NULL, #number of points to sample
                             scale.lambda=c(1,1,1),
                             maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

#записываем в otput
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#находим индексы минимальных значений BIC для каждого значения к
minBICid = apply(BIC,2,which.min)
minBICid
k = 2
min(BIC[,k])

#строим кривую % объясненной разницы
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
devRatMinBIC

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
#Вывод: На основании одного набора омиксных данных вывод о количестве калстеров сделать невозможно



#Bayes
setwd('C:/Users/Irina Makarova/Documents/ADBM/System biology - fri/hw3/fit_bayes_no_thin')
bayfit = tune.iClusterBayes(cpus=1,
                            dt1=mrna_s,
                            dt2=mirna_s,
                            dt3=prot_s,
                            type=c("gaussian","gaussian","gaussian"),
                            K=1:4,
                            n.burnin=1000,
                            n.draw=1200,
                            prior.gamma=c(0.1,0.1,0.1),
                            sdev=0.5)
#сохраняем данные
save.image(file="gbmBayfit.RData")
#load("gbmBayfit.RData")
allBIC = NULL
devratio = NULL
nK = length(bayfit$fit)
for(i in 1:nK){
  allBIC = c(allBIC,bayfit$fit[[i]]$BIC)
  devratio = c(devratio,bayfit$fit[[i]]$dev.ratio)
}
#отображаем кривые deviance ratio и BIC, находим оптимальное количество кластеров
par(mar=c(4.0,4.0,0.5,0.5),mfrow=c(1,2))
plot(1:nK, allBIC,type="b",xlab="k",ylab="BIC",pch=c(1,1,19,1,1,1))
plot(1:nK,devratio,type="b",xlab="k",ylab="Deviance ratio",pch=c(1,1,19,1,1,1))

bayfit$fit[[3]]

best.cluster.Bayes = bayfit$fit[[2]]$clusters
write.table(best.cluster.Bayes, file="clusterMembership.txt",sep='\t',quote=F)

k=2
best.cluster=clusters[,k]
best.fit=output[[k]]$fit[[which.min(BIC[,k])]]
best.fit


#Кластеризация, iClusterPlus, только гены, содержащиеся в С4 cancer modules, n.lambda = NULL
#убираем гены, не содержащиеся в С4 (отбор сделан отдельно в питоне), 85 генов из 200
mrna_sp2 = subset(mrna, 
                  select = -c(ADAMTS4,
                              AMN1,
                              AMPD3,
                              BOC,
                              C11orf52,
                              C18orf1,
                              C1orf162,
                              C1orf38,
                              C4orf34,
                              C6orf192,
                              C7orf55,
                              CAMK2D,
                              CCDC113,
                              CCDC64B,
                              CCDC80,
                              CDCP1,
                              CERCAM,
                              DEPDC6,
                              DNLZ,
                              DOCK8,
                              DTWD2,
                              EIF4EBP3,
                              ELP2,
                              FAM63A,
                              FGD5,
                              FLJ23867,
                              FRMD6,
                              FUT8,
                              FZD4,
                              GIYD2,
                              GK,
                              HDAC11,
                              HIST1H2BK,
                              #HLA-H,
                              HN1,
                              ICAM2,
                              IFIH1,
                              IGSF9,
                              ILDR1,
                              KDM4B,
                              KIF13B,
                              LASS4,
                              LRIG1,
                              LRRC25,
                              MANSC1,
                              MED13L,
                              MEX3A,
                              MRVI1,
                              MTL5,
                              NES,
                              NRARP,
                              OGFRL1,
                              PLCD3,
                              PRKCDBP,
                              PROM2,
                              PSIP1,
                              PVRL4,
                              RBL1,
                              RIN3,
                              RTKN2,
                              S100A16,
                              S1PR3,
                              SAMD9L,
                              SDC4,
                              SEMA5A,
                              SGPP1,
                              SH3KBP1,
                              SHROOM3,
                              SIGIRR,
                              SLC19A2,
                              SLC43A2,
                              SLCO3A1,
                              SNHG1,
                              SNORA8,
                              TIGD5,
                              TRIM45,
                              TTC23,
                              TTC39A,
                              WWC2,
                              YPEL2,
                              ZKSCAN1,
                              ZNF37B,
                              ZNF552,
                              ZNF680,
                              ZNRF3
                  ))


#####корреляция
correlation = cor(mrna_sp2,method = 'spearman')
heatmap(correlation,
        symm = TRUE,
        distfun = function(x){as.dist(1-x)}) 
w = which(abs(correlation)>0.8 & row(correlation)<col(correlation), arr.ind=TRUE)
corr_mrna = matrix(colnames(correlation)[w],ncol=2)
corr_mrna
unique(corr_mrna[,2])
mrna_sp2 = subset(mrna_sp2, select = -c(AKAP12, TCF4, COL6A1, APBB1IP, CCNA2, POSTN,  C1QB, CTSK, JAM3))

####Шкалирование
mrna_ssp2 = apply(mrna_sp2, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

####iCluster
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=1, 
                             dt1=mrna_ssp2,
                             dt2=mirna_ssp,
                             dt3=prot_ssp,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,
                             n.lambda=NULL, #number of points to sample
                             scale.lambda=c(1,1,1),
                             maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

#записываем в otput
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#находим индексы минимальных значений BIC для каждого значения к
minBICid = apply(BIC,2,which.min)
minBICid
k = 2
min(BIC[,k])

#строим кривую % объясненной разницы
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
devRatMinBIC

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
#кривая немного снижается после 3 кластеров, но в целом кажется, что данных не хватает для однозначного вывода. Пробуем использовать другой параметр n.lambda 


#Кластеризация, iClusterPlus, шкалированные данные, часть признаков,только гены, содержащиеся в С6 oncogenic signature gene sets, n.lambda = NULL
#убираем гены, не содержащиеся в С6 (отбор сделан отдельно в питоне), 51 генов из 200
mrna_sp3 = subset(mrna, 
                  select = -c(ADAMTS4,
                              BOC,
                              C18orf1,
                              C1orf162,
                              C1orf38,
                              C4orf34,
                              C6orf192,
                              C7orf55,
                              CCDC113,
                              CCDC64B,
                              CERCAM,
                              DEPDC6,
                              DNLZ,
                              DTWD2,
                              ELP2,
                              FAM63A,
                              FGD5,
                              FLJ23867,
                              GIYD2,
                              HIST1H2BK,
                              HN1,
                              IFITM2,
                              IGSF9,
                              ILDR1,
                              LAPTM4B,
                              LASS4,
                              LYN,
                              MEX3A,
                              MRVI1,
                              MTL5,
                              PLCD3,
                              PLCD4,
                              PREX1,
                              PRKCDBP,
                              PRNP,
                              PROM2,
                              PVRL4,
                              RIN3,
                              SDC1,
                              SEMA4A,
                              SH3KBP1,
                              SIGIRR,
                              SLC43A2,
                              SNHG1,
                              SNORA8,
                              TBC1D4,
                              TIGD5,
                              YPEL2,
                              ZNF37B,
                              ZNF552,
                              ZNF680))


#####корреляция
correlation = cor(mrna_sp3,method = 'spearman')
heatmap(correlation,
        symm = TRUE,
        distfun = function(x){as.dist(1-x)}) 
w = which(abs(correlation)>0.8 & row(correlation)<col(correlation), arr.ind=TRUE)
corr_mrna = matrix(colnames(correlation)[w],ncol=2)
corr_mrna
unique(corr_mrna[,2])
mrna_sp3 = subset(mrna_sp3, select = -c(AKAP12, TCF4, COL6A1, APBB1IP, CCDC80, CCNA2, POSTN, LRRC25, C1QB, CTSK, JAM3))

####Шкалирование
mrna_ssp3 = apply(mrna_sp3, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))

####iCluster
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=1, 
                             dt1=mrna_ssp3,
                             dt2=mirna_ssp,
                             dt3=prot_ssp,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,
                             n.lambda=NULL, #number of points to sample
                             scale.lambda=c(1,1,1),
                             maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

#записываем в otput
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#находим индексы минимальных значений BIC для каждого значения к
minBICid = apply(BIC,2,which.min)
minBICid
k = 2
min(BIC[,k])

#строим кривую % объясненной разницы
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
devRatMinBIC

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
#вывод о количестве кластеров не удается сделать 

clusters=getClusters(output)
rownames(clusters)=rownames(mrna_ssp)
colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
write.table(clusters, file="clusterMembership.txt",sep='\t',quote=F)
clusters




#Ниже эксперименты с n.lambda, которые оказались лишены смысла
#Кластеризация по 3 омиксам, iClusterPlus, шкалированные данные, все признаки, n.lambda = 418
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=1,
                             dt1=mrna_s,
                             dt2=mirna_s,
                             dt3=prot_s,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,
                             n.lambda=418, #number of points to sample, 418 - одна из опций для 3 датасетов, указанная в руководстве
                             scale.lambda=c(1,1,1),
                             maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

#записываем в otput
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#находим индексы минимальных значений BIC для каждого значения к
minBICid = apply(BIC,2,which.min)
minBICid

#строим кривую % объясненной разницы
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
devRatMinBIC

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
#кажется, что модель переобучилась 

#Кластеризация по 3 омиксам, iClusterPlus, шкалированные данные, все признаки, n.lambda = 35
for(k in 1:5){
  cv.fit = tune.iClusterPlus(cpus=1,
                             dt1=mrna_s,
                             dt2=mirna_s,
                             dt3=prot_s,
                             type=c("gaussian","gaussian","gaussian"),
                             K=k,
                             n.lambda=35, #number of points to sample
                             scale.lambda=c(1,1,1),
                             maxiter=20)
  save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
}

#записываем в otput
output=alist()
files=grep("cv.fit",dir())
for(i in 1:length(files)){
  load(dir()[files[i]])
  output[[i]]=cv.fit
}
nLambda = nrow(output[[1]]$lambda)
nK = length(output)
BIC = getBIC(output)
devR = getDevR(output)

#находим индексы минимальных значений BIC для каждого значения к
minBICid = apply(BIC,2,which.min)
minBICid
k = 3
min(BIC[,k]) #-204008.2

#строим кривую % объясненной разницы
devRatMinBIC = rep(NA,nK)
for(i in 1:nK){
  devRatMinBIC[i] = devR[minBICid[i],i]
}
devRatMinBIC

plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
     ylab="%Explained Variation")
#нет ожидаемого спада на кривой, одноначный вывод сделать невозможно.
#Предварительный вывод: значение n.lambda, равное 185, представляется самым подходящим для данных такой размерности.
