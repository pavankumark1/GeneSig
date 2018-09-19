#KDD Project 
#Authors Of code: #Pavan Kumar Komma.


#setting Working directory
setwd("~/GSE17537")


#Sources
source("http://www.bioconductor.org/biocLite.R")

#packages
biocLite("affy")
n
biocLite("limma")
n
biocLite("hgu133plus2cdf")
n
biocLite("hgu133plus2.db")
n

#libraries

library(affy)
library(limma)
library(hgu133plus2cdf)
library(hgu133plus2.db)

#files

fn <- dir() 

dir()

str(fn)

#reading raw data
rawdata<- ReadAffy(filenames = fn, cdfname = "hgu133plus2cdf")

#boxplot of rawdata
boxplot(rawdata, col = "blue")

#normalizing the data 
ndata<- justRMA(filenames = fn, cdfname = "hgu133plus2cdf", verbose = TRUE)

#expression of data
edata <- exprs(ndata)

#boxplot of normalized data
boxplot(edata, col = "red")

#changing the column names of expressed data
colnames(edata) <-  c("NRC","NRC","NRC","NRC","RC","NRC","RC",

        "NRC","NRC","NRC","NRC","NRC","NRC","NRC","NRC","NRC",

        "RC","NRC","NRC","NRC","RC","NRC","RC","RC","NRC","RC",

        "NRC","NRC","NRC","NRC","NRC","RC","RC","NRC","NRC",

        "RC","NRC","NRC","RC","NRC","NRC","NRC","NRC","NRC",

        "RC","NRC","RC","RC","RC","RC","RC","RC","RC","NRC","NRC")
                   
                    
cnames <- colnames(edata) 

#Analysis of DEG's using limma ("differentially expresssed genes")

designmatrix <- model.matrix(~0 + factor(cnames))

colnames(designmatrix) <- c("RC","NRC")

View(designmatrix)

fit <- lmFit(edata, designmatrix)

#Contrast matrix 

conmatrix <- makeContrasts(RC-NRC, levels = designmatrix)

fit <- contrasts.fit(fit,conmatrix)

#View(fit)

#eBayes stataistics application

efit <- eBayes(fit)

siggenes<- topTable(efit,coef = 1, adjust = "fdr", number = length(efit))

DEGS <- subset(siggenes,siggenes$P.Value < 0.05 & abs(siggenes$logFC)> 0.7)


#replacing DEG's rownames with genenames 
degrownames <- rownames(DEGS)

#getting the gene names for the rownames of DEG's  
genenames<- unlist(mget(degrownames,hgu133plus2SYMBOL, ifnotfound = NA))


#binding the geneID's to the DEGS 
#this is done to get some clarity with the genenames and their ID's 
DEGS<- cbind(geneIDs = genenames,DEGS)

#list of DEG's 
View(DEGS)


        #''''''' SVM ANALYSIS''''''''''
        
#subsettting data depending on the DEG's rownames
#for SVM analysis, We will only be using the data present in the "svmdata" matrix
install.packages("e1071")

library(e1071)

svmdata <- subset(edata,rownames(edata) %in% rownames(DEGS))

#rownames(svmdata) <-unlist(mget(rownames(svmdata),hgu133plus2SYMBOL, ifnotfound = NA)) 

#GSE17538 data
setwd("~/CIS 635 project/GSE17538")

getwd()

fn2 <- dir()

str(fn2)

#Reading the rawdata
newGSEdata <- ReadAffy(filenames = fn2, cdfname = "hgu133plus2cdf")

boxplot(newGSEdata, col = "blue")

#expressing and normalizing the GSE17538 data
nGSE <- justRMA(filenames = fn2, cdfname = "hgu133plus2cdf", verbose = TRUE)

eGSE <- exprs(nGSE)

boxplot( eGSE, col = "red")

#changing column names of expressed dataset(eGSE)


for(i in 1: nchar(gsms)) {
  
    if (substr(gsms,i,i) =="0"){ newcnames[i]<-"NRC"}
  
  else if (substr(gsms,i,i) == 1) { newcnames[i] <- "RC" }
  else if (substr(gsms,i,i) == 2) { newcnames[i] <- "NA1"}
  else if (substr(gsms,i,i) == 3) { newcnames[i] <- "NA2"}
  
}
colnames(eGSE) <- newcnames

colnames(eGSE)

eGSE <- eGSE[ , !(colnames(eGSE) %in% c("NA1","NA2"))]

#Identifying DEG's 

svmdata2 <- subset(eGSE,rownames(eGSE) %in% rownames(DEGS))

write.csv(svmdata2, file = "testingset.csv")

#test and training data
trsvmdata <- t(svmdata)
trsvmdata2 <- t(svmdata2)


set.seed(12)

trnr <- row.names(trsvmdata)

trnd <- data.frame(trsvmdata)

y  <- as.factor(trnr)

tesr <- row.names(trsvmdata2)

tesd <- data.frame(trsvmdata2)

x<- as.factor(tesr)

#training sets -- trnr, trnd/y (GSE17537)
#testing sets -- tesr, tesd/x (GSE17538)

svm.error <- matrix(0, nrow = 3, ncol = 3)  ## error will be saved
for(cost in 0:2)     ## cost = 1, 2, 4
{
  for(gamma in (-1):1)   ## gamma = 2^gamma/ncol(x.learn)
  {
    i              <- cost+1   ## index for error matrix
    j              <- gamma+2
    svm.fit        <- svm(trnd, y, cost = 2^cost,    
                          gamma = 1)
    svm.predic     <- predict(svm.fit, newdata = tesd)
    svm.error[i,j] <- mean(svm.predic != x)
  }
}
#getting the predic score and creating a table based on prediction.

predscore <- predict(svm.fit, tesd, type = "decision")
View(predscore)

#expected prediction to the original predication 

tab <- table(pred = svm.predic, true = tesr)
View(tab)


