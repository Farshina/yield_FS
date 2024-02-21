function [avg_rmse_data] = yieldKS_evaluate_features(X,y)
% Set up 5-fold cross-validation
cv = cvpartition(length(y), 'KFold', 10);

% Initialize RMSE vectors
rmse_values_svm = zeros(cv.NumTestSets, 1);
rmse_values_knn = zeros(cv.NumTestSets, 1);
rmse_values_forest = zeros(cv.NumTestSets, 1);

% Perform cross-validation
for fold = 1:cv.NumTestSets
    trainIdx = cv.training(fold);
    testIdx = cv.test(fold);

    X_fold_train = X(trainIdx, :);
    y_fold_train = y(trainIdx);
    X_fold_test = X(testIdx, :);
    y_fold_test = y(testIdx);

    % Train the regression models for each fold
    svmModel = fitrsvm(X_fold_train, y_fold_train,"KernelFunction","rbf","Standardize",true); % SVM
    % knnModel = kNNeighborsRegressor(4,'euclidean','uniform'); % KNN
    knnModel = kNNeighborsRegressor(4,'euclidean','distance'); % KNN
    knnModel = knnModel.fit(X_fold_train, y_fold_train);
    forestModel = fitrensemble(X_fold_train, y_fold_train, 'Method', 'bag'); % RF:: Bootstrap Aggregating

    % Make predictions on the test set for each fold
    y_pred_svm = predict(svmModel, X_fold_test); % SVM
    y_pred_knn = knnModel.predict(X_fold_test); % KNN
    y_pred_knn = reshape(y_pred_knn,[],1);
    y_pred_forest = predict(forestModel, X_fold_test); % RF

    % Calculate RMSE for each fold
    rmse_values_svm(fold) = sqrt(mean((y_fold_test - y_pred_svm).^2));
    rmse_values_knn(fold) = sqrt(mean((y_fold_test - y_pred_knn).^2));
    rmse_values_forest(fold) = sqrt(mean((y_fold_test - y_pred_forest).^2));

end

% Average RMSE over all folds
avg_rmse_svm = mean(rmse_values_svm);
avg_rmse_knn = mean(rmse_values_knn);
avg_rmse_forest = mean(rmse_values_forest);
avg_rmse_data = [avg_rmse_svm,avg_rmse_knn,avg_rmse_forest];
end