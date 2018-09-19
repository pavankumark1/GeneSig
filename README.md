# GeneSig
Study the gene signatures for patients with high and low risk recurrence of colon cancer and predicting their likely outcomes


Two microarray datasets (GSE17537 and GSE17538) were collected from Gene Expression Omnibus database. After preprocessing, both the datasets have been normalized using Linear Models for Microarray Data (LIMMA) method to identify the differentially expressed genes. 109 genes have been selected from both the datasets with genes expressed in GSE17537 dataset being the training data and genes in GSE17538 dataset being the test data. This data is used for Support Vector Machine Analysis for determining the prediction accuracy of recurrent and non-recurrent samples. 
