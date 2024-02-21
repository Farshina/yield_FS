yieldKSdataset_total.m

Input and details: all USDANASS files for yeild data and PRISM files for climate data. Combine these information, get rid of outliers, make dataset. 
Output: climate monthly and annual values for KS-based yield- all districts, irrigation districts, rainfed districts (yielddataset_kansas_monthly_total_1981_2018_0mm.xlsx).


yieldKSfs_mRMR_modelbased.m

Input and details: Output datasets of yieldKSdataset_total.m (yielddataset_kansas_monthly_total_1981_2018_0mm.xlsx). Data is prepared for mRMR, then mRMR performed, best features are selected using 3 models (KNN, SVM, RF), 5 fold CV and corresponding RMSE.
Output: RMSE table of the model performances using features 1,2,3,..n


A separate file for the discussion of the paper: curseofdimentionaity_datatofeatureratio.m