function [RMSEtable] = RMSE_modelbased(X,y,predictorNames,load_filename,dstcode)

% MRMR
[idx,scores] = findmrmr(X,y);

N = length(idx);
n_models = 3; % SVM, KNN, RF
RMSEtable = array2table(zeros(N,n_models+1),'VariableNames',["SVM","KNN","RF","Finalvalue"]);
load_filename
for n = 1:N
    rng("default")
    if mod(n,10) == 0
        n
    end
    tempX = X(:,idx(1:n)); % selected columns only
    [avg_rmse_data] = yieldKS_evaluate_features(tempX,y);
    RMSEtable{n,:} = [avg_rmse_data,mean(avg_rmse_data)];
end

% save error table
originalfilename = "RMSEinfo";
savefilename = get_filenameext(originalfilename,load_filename,dstcode);
save(savefilename + ".mat","RMSEtable","predictorNames","scores","idx")
end
