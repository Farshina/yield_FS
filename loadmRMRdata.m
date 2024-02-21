function [T,target_var,predictorNames,X,y] = loadmRMRdata(load_filename)
% read data from file
Tdata = readtable(load_filename);

% select met. variables only
varnames = {"AgDistrict","County","CountyANSI","Year","AgDistrictCode"};
for i = 1:length(varnames)
    var_idx(i) = find(strcmpi(Tdata.Properties.VariableNames,varnames{i}));
end

% get rid of "unnecessary" variables
T = removevars(Tdata,var_idx);

% define target and predictor variables
target_var = "TONS_ACRE";
varnames = T.Properties.VariableNames;
target_idx = find(strcmpi(varnames,target_var));
predictorNames = varnames;
predictorNames(target_idx) = [];
X = T{:,predictorNames};
y = T{:,target_var};
end
