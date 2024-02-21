clear all
close all
clc

ndatapoints = [3353,252,252,444,356,350,431,347,465,456,2085,1268]; % all-districts, 9 district, irrigation, rainfed
ndatasets = length(ndatapoints);
k = 5; % k-fold crossvalidation
dp_perfold = floor(ndatapoints/k);
n_features = 117;
ratio = zeros(n_features,ndatasets);
for i = 1:n_features
    ratio(i,:) = dp_perfold/i;
end
% plot
figure("Units","normalized","OuterPosition",[0 0 1 1])
plottitle = "Ratio of number of datapoints per fold to number of features used for prediction";
xval = 1:n_features;
lgnd = {"all-districts","district 10","district 20","district 30","district 40","district 50","district 60","district 70","district 80","district 90","irrigation districts","rainfed districts"};
markers = '*.........xx';
for j = 1:ndatasets
    semilogy(xval,ratio(:,j),"Marker",markers(j))
    hold on
end
yline(1,"LineStyle","--","Color",'r')
xlim([min(xval) max(xval)])
xlabel("No. of features")
ylabel("Ratio = datapoints per fold ÷ No. of features")
legend(lgnd)
grid minor
title(plottitle)

% save
saveas(gcf,[plottitle + ".png"])
savefig(plottitle)
