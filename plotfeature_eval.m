function [best_n_features] = plotfeature_eval(RMSEtable,load_filename,dstcode)
n_features = height(RMSEtable);
xval = 1:n_features;
n_models = width(RMSEtable);
lgnd = RMSEtable.Properties.VariableNames;
lgnd{end} = "Average Value";
plottitle_temp = "Evaluation Score (RMSE) for Different Prediction Models";
plottitle = get_filenameext(plottitle_temp,load_filename,dstcode);

% moving average
dim = 1; % column
kb = 3;
kf = 3;
RMSE_moving_average = movmean(RMSEtable{:,:}, [kb kf], dim, "omitnan","Endpoints","fill");
best_n_features = 0;

% Plot evaluation scores
figure('units', 'inches', 'outerposition', [0 0 10 6])
markers = '*x+o';
colors = {"#009DD1";"#8a64d6";"#269C7D";"#FF0000"};
for i = 1:n_models
    line_width = 1;
    line_style = "--";
    if i == n_models
        line_width = 1;
        line_style = "-";
    end
    plot(xval,RMSE_moving_average(:,i),"Color",colors{i},"Marker",string(markers(i)),"LineStyle",line_style,"LineWidth",line_width)
    hold on
end
for extra = n_models+1:1:length(colors)
    yline(rand(1,1),"Color",colors{extra},'LineWidth',3)
    lgnd{end+1} = colors{extra};
end
% find the first local minima
TF = islocalmin(RMSE_moving_average(:,end));
%%%%% plot all local minima in the figure
% plot(xval(TF),A(TF),'r*',"MarkerSize",12)
% lgnd{end+1} = "localminima";
% get the first local minima
best_n_features = find(TF == 1, 1, 'first');
bestval = RMSE_moving_average(best_n_features,end);
xline(best_n_features,'--',{"best no. of predictors = " + string(best_n_features), "best RMSE = " + string(round(bestval,2))},'LabelVerticalAlignment','bottom','LabelOrientation','horizontal')
xlabel("No. of Predictors")
ylabel("RMSE")
ylim([0 1.5])
legend(lgnd,"Location","southeast")
title(plottitle)
grid minor

% save figure
originalfilename = "modelRMSEscore";
savefilename = get_filenameext(originalfilename,load_filename,dstcode);
saveas(gcf,[savefilename + ".png"])
savefig(savefilename)
end
