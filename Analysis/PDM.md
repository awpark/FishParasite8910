---
title: "Parasite species distribution modeling"
author: "Tad Dallas, Andrew Park, and John M. Drake"
output:
  html_document:
    fig_height: 6
    fig_width: 6
    highlight: tango
    theme: journal
  pdf_document: default
bibliography: carp.bib
---





## Introduction

### Species distirbution models and the importance of them



### Knowledge gap

What constrains the range of hosts that a parasite can infect? Is there a simple range of host functional traits that can determine the likelihood that a parasite infects a given host species? How well can we predict parasite occurrences given _only_ some host life history traits?

### Thesis paragraph

Here, I apply a series of predictive models in order to predict parasite occurrence across a range of potential host species for a large set of parasites of freshwater fish, using host functional traits, and geographic location. 




### thorns in my side:

* Absence data aren't true absences. Should I even train on these data if the model treats them as true absences?

* How much time to invest reading density estimation literature? 




## Methods

### Data and processing
 We use an existing global database of fish-parasite associations [@strona2013] consisting of over 38000 parasite records spanning a large diversity of parasites (Acanthocephala, Cestoda, Monogenea, Nematoda, Trematoda). In order to allow for cross-validation and accurate prediction, we constrained our ananlyses to parasites with a minimum of 50 host records. In other words, we only examined parasites that had been recorded more than 50 times, but these occurrences could be on fewer than 50 host species. The inclusion of duplicate occurrences was only permitted if the parasite was recorded on a host in a different location, based on latitude and longitude values. Our response variable was parasite occurrence (binary), and was predicted using only host life history traits, and geographic location of host capture. Host trait information was obtained through the FishPest database [@strona2012; @strona2013], and FishBase [@froese2010]. Host traits descriptions are provided in Table 1
 
 
### Model formulation 
  We trained a series of models in order to compare predictive performance of different techniques. Each model was trained on 70% of the data, and accuracy was determined from the remaining 30%. This process was repeated $z$ times ($z$ = 20). We generated background data by randomly sampling host species where parasite $i$ was not recorded. To maintain proportional training data, the number of random samples was selected to be five times greater than the occurrence records. 
  
### Models used
 discuss null predictions scenario, and then go into other algorithms used (brt, svm, lr, rf)







### Strona's FishPest database


```r
load('/media/drakelab/Lexar/8910Project/Analysis/pest.Rdata')
library(rfishbase); library(gbm); library(ROCR); library(e1071); library(randomForest)
## Creating a `fishMatrix` equiv.
#edgePest=edge2matrix(cbind(pest[,'H_SP'], pest[,'P_SP']))
#colnames(edgePest)=unique(pest[,'P_SP'])
#rownames(edgePest)=unique(pest[,'H_SP'])
#edgePest=edgePest[,-which(colSums(edgePest) < 50)]

# storage for model outputs
baseline.auc=vector(); 
brtModel=list(); brt.best.iter=list(); brt.preds=list(); brt.perf=list(); brt.perfAUC=vector()
svmModel=list(); svm.preds=list(); svm.perf=list(); svm.perfAUC=vector()
lrModel=list(); lr.preds=list(); lr.perf=list(); lr.perfAUC=vector()
rfModel=list(); rf.preds=list(); rf.perf=list(); rf.perfAUC=vector()


## Creating a pointsObject
#hostPest = unique(pest[,-c(1:3)], MARGIN=1)
hostPest = pest[,-c(1:3)]
parPest = names(summary(pest[,'P_SP'])>50); parPest=parPest[-100]

ret=list()
reps=20
for(z in 1:reps){
  for(i in 1:length(parPest)){
#create presence vector
  presence=rep(0, nrow(hostPest))
  presence[which(hostPest[,'P_SP'] == parPest[i])]=1

#make some 'na' into actual NAs
  ugh=which(hostPest=='na', arr.ind=TRUE)
  hostPest[ugh]=NA
  hostPest[,7:17]=apply(hostPest[,7:17],2, as.numeric)

#Only train on some of the absences, since they're totally not true absences
  cutDown=c(which(presence==1), sample(which(presence==0), 5*sum(presence)))
  presence1=presence[cutDown]
  dat=hostPest[cutDown, -c(1,5)]

#Impute the data
  impDat=rfImpute(dat[,-c(1:3)], presence1)
  dat=impDat[,-1]

flag=0
 while(flag == 0){
  # This makes sure that the test set contains at least 4 hosts on which the parasite actually occurs
  inds=sample(1:nrow(dat), 0.7*nrow(dat))
  if(sum(presence1[inds]) < 4){inds[1:4] = which(presence1 == 1)[1:4]}

  #Set up a prelim train set and a test set
  train = dat[inds,]
  test = dat[-inds,]  
  if(all(unique(train$GEO) %in% unique(test$GEO))){flag=1}
 }
 
 #Presences
  prezTR=presence1[inds]
  prezTE=presence1[-inds]

  
  #baseline expectations and null models
  baseline.auc[i] = performance(prediction(sample(presence[-inds], length(presence[-inds])), presence[-inds]), 'auc')@y.values
   
 ##trained models
  #boosted regression trees
    brtModel[[i]] = gbm(prezTR ~ ., data=train, n.trees=50000, interaction.depth=4, distribution='bernoulli')
    brt.best.iter[[i]] = gbm.perf(brtModel[[i]], method="OOB")
    brt.preds[[i]] = prediction(predict(brtModel[[i]], newdata=test, n.trees=brt.best.iter[[i]]), prezTE)
    brt.perf[[i]] = performance(brt.preds[[i]],"tpr","fpr")
    brt.perfAUC[i]=unlist(performance(brt.preds[[i]], 'auc')@y.values)

  #support vector machines
    svmModel[[i]] = svm(prezTR ~ ., data=train, probability=TRUE)
    svm.preds[[i]] = prediction(predict(svmModel[[i]], test), prezTE)
    svm.perf[[i]] = performance(svm.preds[[i]],"tpr","fpr")
    svm.perfAUC[i] = unlist(performance(svm.preds[[i]], 'auc')@y.values)

  #logistic regression
    lrModel[[i]] = glm(prezTR ~ ., data=train, family=binomial)
    lr.preds[[i]] = prediction(predict(lrModel[[i]], test), prezTE)
    lr.perf[[i]] = performance(lr.preds[[i]],"tpr","fpr")
    lr.perfAUC[i] = unlist(performance(lr.preds[[i]], 'auc')@y.values)

  #random forest
    rfModel[[i]] = randomForest(prezTR ~ ., data=train)
    rf.preds[[i]] = prediction(predict(rfModel[[i]], test), prezTE)
    rf.perf[[i]] = performance(rf.preds[[i]],"tpr","fpr")
    rf.perfAUC[i] = unlist(performance(rf.preds[[i]], 'auc')@y.values)
  print(i)
  }

 ret[[z]]=cbind(baseline.auc, brt.perfAUC, svm.perfAUC, lr.perfAUC, rf.perfAUC)
 colnames(ret[[z]])=c('BASE', 'BRT', 'SVM', 'LR', 'RF')
 }
```






## Results









## Discussion



## Tables

Table 1: Host traits examined and their corresponding units

Table 2: Details for parasites modeled, including number of occurrences, and life history information.


## Figures








## References 




















