clear all
close all
clc

% load nbest_features data
load("best_n_features_all_dataset.mat")

% all districts dataset
load_filename = "yielddataset_kansas_monthly_total_1981_2018_0mm.xlsx";
best_n_features_idx = 1; % init index
dstcode = NaN;

% prepare data for mRMR
[T,target_var,predictorNames,X,y] = loadmRMRdata(load_filename);

% perform mRMR
[idx,scores] = findmrmr(X,y);

% select features
best_n_features = best_n_features_array(best_n_features_idx);
[selectedvars,selected_scores] =  select_features_featurescores(best_n_features,scores,idx,predictorNames);
best_n_features_idx = best_n_features_idx + 1; % increment index


% plot selected features
plotmRMR_selectedvars(selected_scores, selectedvars, load_filename,dstcode);

% district-specific dataset
Tdata = readtable(load_filename);
dstcode_all = unique(Tdata.AgDistrictCode);
for i = 1:length(dstcode_all)
    dstcode = dstcode_all(i);
    dst_idx = find(Tdata.AgDistrictCode == dstcode);

    % extract district-specific table
    T_dst = T(dst_idx,:);
    X_dst = T_dst{:,predictorNames};
    y_dst = T_dst{:,target_var};

    % perform mRMR
    [idx,scores] = findmrmr(X_dst,y_dst);

    % select features
    best_n_features = best_n_features_array(best_n_features_idx);
    [selectedvars,selected_scores] =  select_features_featurescores(best_n_features,scores,idx,predictorNames);
    best_n_features_idx = best_n_features_idx + 1; % increment index


    % plot selected features
    plotmRMR_selectedvars(selected_scores, selectedvars, load_filename,dstcode);
end

% set dstcode to nan
dstcode = NaN;

% irrigation districts dataset
load_filename = "yielddataset_irrigationdst_kansas_monthly_total_1981_2018_0mm.xlsx";

% prepare data for mRMR
[T,target_var,predictorNames,X,y] = loadmRMRdata(load_filename);

% perform mRMR
[idx,scores] = findmrmr(X,y);

% select features
best_n_features = best_n_features_array(best_n_features_idx);
[selectedvars,selected_scores] =  select_features_featurescores(best_n_features,scores,idx,predictorNames);
best_n_features_idx = best_n_features_idx + 1; % increment index


% plot selected features
plotmRMR_selectedvars(selected_scores, selectedvars, load_filename,dstcode);

% rainfed districts dataset
load_filename = "yielddataset_rainfeddst_kansas_monthly_total_1981_2018_0mm.xlsx";

% prepare data for mRMR
[T,target_var,predictorNames,X,y] = loadmRMRdata(load_filename);

% perform mRMR
[idx,scores] = findmrmr(X,y);

% select features
best_n_features = best_n_features_array(best_n_features_idx);
[selectedvars,selected_scores] =  select_features_featurescores(best_n_features,scores,idx,predictorNames);
best_n_features_idx = best_n_features_idx + 1; % increment index


% plot selected features
plotmRMR_selectedvars(selected_scores, selectedvars, load_filename,dstcode);

function [selectedvars,selected_scores] =  select_features_featurescores(best_n_features,scores,idx,predictorNames)
% select best n features and their respective scores
selectedvars = predictorNames(idx(1:best_n_features));
selected_scores = scores(idx(1:best_n_features));

% add a very small number to the very small scores (just for the plotting purpose)
% selected score is alrady sorted!!!
for i = 2:length(selected_scores)
    if selected_scores(i-1) - selected_scores(i) > 0.01
        selected_scores(i) = selected_scores(i) + 0.01;
    end
end
end


function plotmRMR_selectedvars(selected_scores, selectedvars, load_filename,dstcode)
originalfilename = "MRMR KS (selected features)";
savefilename = get_filenameext(originalfilename,load_filename,dstcode);
plottitle = savefilename;

% Plot scores
figure('units', 'inches', 'outerposition', [0 0 10 6])
bar(selected_scores)
xlabel("Predictor rank")
ylabel("Predictor importance score (Normalized)")
xticks(1:length(selectedvars))
xticklabels(strrep(selectedvars,"_","\_"))
xtickangle(90)
title(plottitle)

% saveplot
saveas(gcf,[savefilename + ".png"])
savefig(savefilename)
% close
end
