#At the top of script needed libraries are imported 
#then the data is scaled and centered  before being split into train and test
#After being spilt into train and test , differential expression is performed.
#p-values from differential expression are shown in a histogram.
#Hits from differential expression are then fed into a random forest classfier and examined for variable importance.
#The top 20 genes are then mapped from probes to gene symboles and about output for viewing and futher work.
#then the 20  genes are used to prepare a Random Forest classfier  and model performance is reported in terms of ROC and confusion matrix.
# note to examine another dataset place the object contianing your dataset in the data function on line 15.
