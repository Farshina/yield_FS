clear all
close all
clc
%%
% all districts dataset
load_filename = "yielddataset_kansas_monthly_total_1981_2018_0mm.xlsx";
dstcode = NaN;

% % prepare data for mRMR
% [T,target_var,predictorNames,X,y] = loadmRMRdata(load_filename);
% height(T)
% % perform mRMR and evaluate RMSE with 4 models
% RMSEtable_all = RMSE_modelbased(X,y,predictorNames,load_filename,dstcode);

% % % district-specific dataset
% Tdata = readtable(load_filename);
% dstcode_all = unique(Tdata.AgDistrictCode);
% for i = 1:length(dstcode_all)
%     dstcode = dstcode_all(i);
%     dst_idx = find(Tdata.AgDistrictCode == dstcode);
% 
%     % extract district-specific table
%     T_dst = T(dst_idx,:);
%     X_dst = T_dst{:,predictorNames};
%     y_dst = T_dst{:,target_var};
%     height(T_dst)
%     % % Perform mRMR and evaluate RMSE with 4 models
%     % RMSEtable_district = RMSE_modelbased(X_dst,y_dst,predictorNames,load_filename,dstcode);
% end
% 
% % set dstcode to nan
% dstcode = NaN;
% 
% % irrigation districts dataset
% load_filename = "yielddataset_irrigationdst_kansas_monthly_total_1981_2018_0mm.xlsx";
% 
% % prepare data for mRMR
% [T,target_var,predictorNames,X,y] = loadmRMRdata(load_filename);
% height(T)
% 
% % % perform mRMR and evaluate RMSE with 4 models
% % RMSEtable_irr = RMSE_modelbased(X,y,predictorNames,load_filename,dstcode);
% 
% % rainfed districts dataset
% load_filename = "yielddataset_rainfeddst_kansas_monthly_total_1981_2018_0mm.xlsx";
% 
% % prepare data for mRMR
% [T,target_var,predictorNames,X,y] = loadmRMRdata(load_filename);
% height(T)

% % perform mRMR and evaluate RMSE with 4 models
% RMSEtable_rain = RMSE_modelbased(X,y,predictorNames,load_filename,dstcode);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% shortcut: load all RMSEtable files
% files = dir(".\RMSEinfo*.mat");
% best_n_features_array = [];
% for i = 1:length(files)
%     load_filename = files(i).name;
%     load(load_filename)
%     % get dst
%     pat =" district " + digitsPattern(1) + "0";
%     dstcode = NaN;
%     if contains(load_filename,pat)
%         dstcode = extractBetween(load_filename,"district ",".mat");
%         dstcode = str2double(dstcode);
%     end
%     % plot RMSE scores
%     [best_n_features] = plotfeature_eval(RMSEtable,load_filename,dstcode);
%     best_n_features_array = [best_n_features_array; best_n_features];
% end
% 
% % savedata
% save("best_n_features_all_dataset.mat","best_n_features_array")
% 
% 
