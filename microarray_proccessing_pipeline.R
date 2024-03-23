#importing libraries 
library(cancerclass)
library(golubEsets)
library(dplyr)
library(cluster)
library(hu6800.db)
library(tidyverse)
library(limma)

library(randomForest)
library(varImp)

#reading in and organizing data

data(Golub_Merge)
df<-data.frame(Golub_Merge$T.B.cell)
df2<-data.frame(Golub_Merge$Treatment)
df3<-data.frame(Golub_Merge$ALL.AML)
Main_df<-data.frame(exprs(Golub_Merge))
Main_df<-t(Main_df)
Main_df<-data.frame((Main_df))

Main_df$cancer_type<-df3$Golub_Merge.ALL.AML
library(caret)
set.seed(1)

#use 70% of dataset as training set and 30% as test set
training_rows<-createDataPartition(y=Main_df$cancer_type,p = .7,list = FALSE)

#sample <- sample(c(TRUE, FALSE), nrow(Main_df_TOP20), replace=TRUE, prob=c(0.7,0.3))
train_AA  <- Main_df[training_rows, ]
test   <- Main_df[-training_rows, ]
Genes<-train_AA %>% dplyr::select(-c('cancer_type'))
#Data_analysis part 1 diffenerial expression 


library("limma")

fit <- lmFit(t(Genes), design=model.matrix(~train_AA$cancer_type))

fit <- eBayes(fit)
tt <- topTable(fit, coef=2,number = ncol(Genes))
tt$names<-row.names(tt)
adjusted_P.values<-tt$adj.P.Val

#plot of the adjusted pvalues, many genes appear to be differentialy expressed

hist(adjusted_P.values,col = "green",main = "Differential Expression Analysis Histogram",xlab = "Adjusted P-Value",ylab = "number of genes")

#isolating differential  expressed genes ,making a table for viewing.

tt_sig<-tt[tt$adj.P.Val <= .05,]

write.csv(tt_sig,"C:/Users/ccape/OneDrive/Documents/difgenes.csv")



#recoupling differntialaly  expressed genes to the the phenotype (ALL or AML)
Diff_genes<-Genes[,(row.names(tt_sig))]
train_AA<-t(train_AA)
train_AA<-data.frame(train_AA)

train_AA<-train_AA[c(names(Diff_genes),"cancer_type"),]
train_AA<-t(train_AA)
train_AA<-data.frame(train_AA)

#data analysis part 2 sorting the differential expressed genes by variable importance  via random forest.


library(randomForest)
library(varImp)


model <- randomForest( as.factor(train_AA$cancer_type)~.,data=train_AA,importance=TRUE)
importance_df<-data.frame(importance(model))
importance_df$genenames<-row.names(importance_df)
importance_df<-importance_df[order(importance_df$ALL, decreasing = TRUE),]

#visualizing the spread of importance of predicting ALL via a histogram

hist(importance_df$ALL,col = "blue",main = "Feature Importance Histogram",xlab = "feature importance",ylab = "number of genes")
#repeating for the spread of importance of predicting AML.
hist(importance_df$AML,col="violet")
#As shown in the plots many of the differenial expressed genes are actually not useful in predicting phenotype. 
#this may be due to collinearity. Therefore we will narrow down our hits to the top 50 genes for prediciting ALL since the authors did so.


Top20_RFgenes<-importance_df$genenames[1:20]
train_AA<-t(train_AA)
train_AA<-data.frame(train_AA)
train_AA<-train_AA[c(Top20_RFgenes,"cancer_type"),]
train_AA<-t(train_AA)
train_AA<-data.frame(train_AA)


Top_Genes<-train_AA %>% dplyr::select(-'cancer_type')
Top_Genes_mat<-matrix(as.numeric(as.matrix(Top_Genes)),nrow = nrow(Top_Genes),ncol =ncol(Top_Genes) )
colnames(Top_Genes_mat)<-colnames(Main_df_TOP20 %>% dplyr::select(-'cancer_type'))
Top_genes_mat<-scale(Top_Genes_mat)
New_df<-as.data.frame(Top_Genes_mat)
New_df$cancer_type<-train_AA$cancer_type

#annotation and outputing  a CSV to aid in text mining/literature searching.

library("AnnotationDbi")


PROBES<- as.character(Top20_RFgenes)
OUT_ <- na.omit(AnnotationDbi::select(hu6800.db,keys= PROBES, columns=c("SYMBOL", "ENTREZID", "GENENAME"),keytype="PROBEID"))

write.csv(OUT_,"C:/Users/ccape/OneDrive/Documents/Genes_rev.csv")

#removing genes not supported by lit


supported_genes<-read.csv("C:/Users/ccape/OneDrive/Documents/supported_Genes.csv")

#writing a text file for annoation 
write.table(supported_genes$SYMBOL,"C:/Users/ccape/OneDrive/Documents/supported_gene_names.txt")





library(caret)


RF_model<-randomForest(as.factor(New_df$cancer_type)~.,data =New_df,importance=TRUE)

prediction <- predict(RF_model, test, 
                      type="response")


confuse <- confusionMatrix(data=as.factor(prediction), reference = as.factor(test$cancer_type))
con_matrix<-data.frame(confuse$byClass)

write.csv(con_matrix,"C:/Users/ccape/OneDrive/Documents/result20genes.csv")




library (pROC)

roc_object <- roc( as.factor(test$cancer_type),as.numeric(as.factor(prediction)))

plot(roc_object)
auc(roc_object)

set.seed(1)



Top22_RFgenes<-importance_df$genenames[1:22]
train_AA<-t(train_AA)
train_AA<-data.frame(train_AA)
train_AA<-train_AA[c(Top20_RFgenes,"cancer_type"),]
train_AA<-t(train_AA)
train_AA<-data.frame(train_AA)


Top_Genes<-train_AA %>% dplyr::select(-'cancer_type')
Top_Genes_mat<-matrix(as.numeric(as.matrix(Top_Genes)),nrow = nrow(Top_Genes),ncol =ncol(Top_Genes) )
colnames(Top_Genes_mat)<-colnames(Main_df_TOP20 %>% dplyr::select(-'cancer_type'))
Top_genes_mat<-scale(Top_Genes_mat)
New_df<-as.data.frame(Top_Genes_mat)
New_df$cancer_type<-train_AA$cancer_type


PROBES<- as.character(Top22_RFgenes)
OUT_rev<- na.omit(AnnotationDbi::select(hu6800.db,keys= PROBES, columns=c("SYMBOL", "ENTREZID", "GENENAME"),keytype="PROBEID"))

write.csv(OUT_rev,"C:/Users/ccape/OneDrive/Documents/Genes_rev2.csv")





RF_model2<-randomForest(as.factor(New_df$cancer_type)~.,data =New_df,importance=TRUE)

prediction2 <- predict(RF_model, test, 
                      type="response")


confuse2 <- confusionMatrix(data=as.factor(prediction), reference = as.factor(test$cancer_type))
con_matrix2<-data.frame(confuse2$byClass)

write.csv(con_matrix,"C:/Users/ccape/OneDrive/Documents/result20genes.csv")




library (pROC)

roc_object <- roc( as.factor(test$cancer_type),as.numeric(as.factor(prediction)))

plot(roc_object)
auc(roc_object)




