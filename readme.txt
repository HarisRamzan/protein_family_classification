1:  pdb_data_no_dups and pdb_data_seq are files which contains sequences,structure information and labels.
2:  both files are merged after cleaning data and by considering structureid columns in both files by taking join and named as FinalDS.csv
3:  After removing redundant sequences and labels cleaning (hydrolase,transferese,other) file created as PreprocessedDataset.csv
4:  After extracting features files is saved as FeatureExtracted.csv
5:  2 files named as svm and xgboost contains the code to train and test the model and evaluation matric f1 score and accuracy are calculated.


*******Note******
>>>>File name as Featureextraction takes about 10 minutes to clean,preprocess, and cleanzing the data according to the requirement.
>>>> To save time only run xgboost and svm named files by following the proper file directory.



#dataset can be found at https://www.kaggle.com/shahir/protein-data-set
