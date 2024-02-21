clear all
close all
clc

%% Load USDA-NASS yield data
path = "D:\Academics\Project_PRECISION_AG_big_files\yield-KS\data and code\USDANASS_1981_2018";
files = dir("USDANASS_1981_2018\*.csv");

% USDA-NASS yielddata 1981-2018
[T,vars] = load_USDANASSdata(path, files);

% create main table
T_main = table(T.Year,VariableNames="Year"); %%%%%%%%%%%%%%%%%% T_main is the main table

% only keep the "useful" columns with unique values
[T, T_main] = clean_table(T, T_main, vars);

clear T vars

% group counties as per district code
[T_county, cnt_unique] = groupcounties(T_main); %%%%%%%%%%%%%%%%%% T_county is the table for each county

% From table, remove for "other" counties
idx_other_tmain = find(strcmpi(T_main.County,"OTHER (COMBINED) COUNTIES"));
T_main(idx_other_tmain,:) = [];

% define table type ['total' or 'average']
tabletype = 'Total';
%% Load PRISM_data into table

path = "D:\Academics\Project_PRECISION_AG_big_files\yield-KS\data and code\PRISM_1981_2018";
files = dir("PRISM_1981_2018\*.csv");

% definne column names
PRISM_keywords = {'ppt'	'tmin'	'tmean'	'tmax'	'tdmean'	'vpdmin'	'vpdmax'};

% iterate through each county
for i = 1:length(files)
    filename = files(i).name;
    fullfilename = fullfile(path,filename);
    currentcounty = extractBefore(filename,".csv");

    % all other parameters
    t = readtable(fullfilename);
    t_years = year(t.Date); % extract years for that county
    t_yr_unique = unique(t_years);
    t_vars = t.Properties.VariableNames;

    % initialize the variable names
    currentyear = t_yr_unique(1);
    yr_idx = find(t_years == currentyear);
    PRISM_yr_county_date = t{yr_idx,t_vars{1}}; % extract the date idx
    PRISM_month_idx = month(PRISM_yr_county_date); % extract the month idx
    month_names = unique(month(PRISM_yr_county_date,'shortname'),'stable');
    if i == 1
        T_main = PRISM_monthly_varnames(month_names,T_main,PRISM_keywords,tabletype);
    end

    % iterate through each year
    for j = 1:length(t_yr_unique)
        currentyear = t_yr_unique(j);
        yr_idx = find(t_years == currentyear);

        % PRISM year/county-specific data
        PRISM_yr_county_tdata = zeros(length(yr_idx),length(t_vars) - 1);

        for t_data_var_idx = 1:length(t_vars) - 1
            PRISM_yr_county_tdata(:,t_data_var_idx) = t{yr_idx,t_vars{t_data_var_idx + 1}}; % t_data_idx + 1 cause we want to avoid the date
        end

        PRISM_yr_county_date = t{yr_idx,t_vars{1}}; % extract the date idx
        PRISM_month_idx = month(PRISM_yr_county_date); % extract the month idx
        month_names = unique(month(PRISM_yr_county_date,'shortname'),'stable');

        % find which row in T_main corresponds that county that year
        [T_main_county_yr_idx] = find_T_main_row_idx(T_main,currentcounty,currentyear);

        if ~isempty(T_main_county_yr_idx) % if the corresponding row exists in T_main
            % assign yearly or monthly total/avg data
            [T_main] = T_main_assign_data(T_main,PRISM_keywords,T_main_county_yr_idx,PRISM_yr_county_tdata,month_names,PRISM_month_idx);
            % insert special variables for no ppt days- both annual and monthly
            [T_main] = T_main_assign_nopptdays(T_main,tabletype,month_names,T_main_county_yr_idx,PRISM_yr_county_tdata,PRISM_month_idx,i,j);
            % insert special variables for growing degree days- both annual and monthly
            [T_main] = T_main_assign_GDD(T_main,tabletype,month_names,T_main_county_yr_idx,PRISM_yr_county_tdata,PRISM_month_idx,i,j);
        end

    end
end

%% split irrigation (10~60) and rainfed (70~90) districts
% irrigation 10-60; rainfed 70-90
irr_dst = 10:10:60;
irrigation_idx = [];
for i = 1:length(irr_dst)
    irrigation_idx = [irrigation_idx; find(T_main.AgDistrictCode == irr_dst(i))];
end
T_irrigation = T_main(irrigation_idx,:);

rainfed_dst = [70,80,90];
rainfed_idx = [];
for i = 1:length(rainfed_dst)
    rainfed_idx = [rainfed_idx; find(T_main.AgDistrictCode == rainfed_dst(i))];
end
T_rainfed = T_main(rainfed_idx,:);
%% Save file
filename = 'yielddataset_kansas_monthly_total_1981_2018_0mm.xlsx';
writetable(T_main,filename);
filename = 'yielddataset_irrigationdst_kansas_monthly_total_1981_2018_0mm.xlsx';
writetable(T_irrigation,filename);
filename = 'yielddataset_rainfeddst_kansas_monthly_total_1981_2018_0mm.xlsx';
writetable(T_rainfed,filename);
%% All functions
function [T,vars] = load_USDANASSdata(path, files)
C = {};
for i = 1:length(files)
    filename = files(i).name;
    fullfilename = fullfile(path,filename);
    t = readtable(fullfilename);
    C{end+1} = t;
end

T = vertcat(C{:});

% rename value to tons/acre
vars = T.Properties.VariableNames;
idx = strncmpi(vars,'Value',2);
vars(idx) = {'TONS / ACRE'};
end

function [T, T_main] = clean_table(T, T_main, vars)
% find the unique cols and get rid of NaN / columns with the same value for
% all input. only keep the "useful" columns with unique values
for i = 1:length(vars)
    if isnumeric(T{:,i})
        T{:,i}(isnan(T{:,i})) = -99999; %%%%%%%%%%%%%%%%%% each nan is considered as a unique value (i know right? :/)
    end
    unique_val{i} = unique(T{:,i}); %%%%%%%%%%%%%%%%%% see if the col has all same value (i.e., "Kansas")
    if length(unique_val{i}) > 1
        name = vars{i};

        T_main.(name) = T{:,i};
    end
end
end

function [T_county, cnt_unique] = groupcounties(T_main)
% group counties as per district code
% because of "others (counties)" in USDA-NASS files, this step is needed
countynames = T_main.County;
cnt_unique = unique(countynames);
dstcodes = zeros(length(cnt_unique),1);
T_county = table(cnt_unique,dstcodes,VariableNames={'County','AgDistrictCode'});

for i = 1:length(cnt_unique)
    currentcounty = string(cnt_unique(i));
    currentidx = find(T_main.County == currentcounty);
    currentdstcodes = unique(T_main.AgDistrictCode(currentidx));
    if length(currentdstcodes) == 1
        T_county.AgDistrictCode(i) = currentdstcodes;
    end
end
end

function [T_main] = PRISM_monthly_varnames(month_names,T_main,PRISM_keywords,tabletype)
% initialize the variable names and reserve a column for them in T_main
nrows = length(T_main.Year);
for PRISM_var = 1:length(PRISM_keywords)
    T_main_varname = horzcat(PRISM_keywords{PRISM_var},'_annual_',tabletype);
    T_main{:,T_main_varname} = zeros(nrows,1);
    for current_month = 1:length(month_names)
        T_main_varname = horzcat(PRISM_keywords{PRISM_var},'_',month_names{current_month},'_',tabletype);
        T_main{:,T_main_varname} = zeros(nrows,1);
    end
end
end

function [T_main_county_yr_idx] = find_T_main_row_idx(T_main,currentcounty,currentyear)
% find which row in T_main corresponds that county that year
T_main_county_idx = find(strcmpi(T_main.County,currentcounty));
T_main_yr_idx = find(T_main.Year == currentyear);
T_main_county_yr_idx = intersect(T_main_county_idx,T_main_yr_idx); % here to be noted that in T_main all year does not exist for all county
if isempty(T_main_county_yr_idx)
    text_out = "County: " + currentcounty + ", Year: " + string(currentyear) + " yield info does not exist in T_main. \n";
else
    text_out = "County: " + currentcounty + ", Year: " + string(currentyear) + " yield info located in T_main. Row: " + string(T_main_county_yr_idx) + "\n";
end
% append county year yieldinfo in text file
fid = fopen('county_year_yieldinfo.txt', 'a+');
fprintf(fid, text_out);
fclose(fid);
% fprintf(text_out)
end

function [T_main] = T_main_assign_data(T_main,PRISM_keywords,T_main_county_yr_idx,PRISM_yr_county_tdata,month_names,PRISM_month_idx)
% iterate thorough all column
T_main_varnames = T_main.Properties.VariableNames;
var_startidx = find(strcmpi(T_main_varnames, 'TONS / ACRE')) + 1; % start index
var_endidx = length(T_main_varnames);

for var_idx = var_startidx:var_endidx
    current_var = T_main_varnames{var_idx};
    for PRISM_key_idx = 1:length(PRISM_keywords)
        if contains(current_var, PRISM_keywords(PRISM_key_idx))
            break % all we need here is the PRISM_key_idx
        end
    end
    if contains(current_var, "annual") % if the data is annual
        T_main{T_main_county_yr_idx,current_var} = sum(PRISM_yr_county_tdata(:,PRISM_key_idx)); % total annual data for currentcounty currentyear
    else % if the data is monthly
        for current_month = 1:length(month_names)
            if contains(current_var, month_names(current_month))
                break % all we need here is the current_month_idx
            end
        end
        month_idx = find(PRISM_month_idx == current_month);
        T_main{T_main_county_yr_idx,current_var} = sum(PRISM_yr_county_tdata(month_idx,PRISM_key_idx)); % total monthly data for currentcounty currentyear
    end
end
end

function [T_main] = T_main_assign_nopptdays(T_main,tabletype,month_names,T_main_county_yr_idx,PRISM_yr_county_tdata,PRISM_month_idx,i,j)
% insert special variables for no ppt days- both annual and monthly
nrows = length(T_main.Year);
if i == 1 && j == 1
    T_main_varname = horzcat('nopptdays_annual_',tabletype);
    T_main{:,T_main_varname} = zeros(nrows,1);
    for current_month = 1:length(month_names)
        T_main_varname = horzcat('nopptdays_',month_names{current_month},'_',tabletype);
        T_main{:,T_main_varname} = zeros(nrows,1);
    end
end

T_main_varnames = T_main.Properties.VariableNames;
var_startidx = find(contains(T_main_varnames, 'nopptdays_annual_')); % start index
var_endidx = length(T_main_varnames);

for var_idx = var_startidx:var_endidx
    current_var = T_main_varnames{var_idx};
    PRISM_key = 'ppt';
    PRISM_key_idx = 1; % hard-coded:::::::: index of ppt is 1
    % ppt_thr = 0.196850394; % 5 mm = 0.196850394 inches
    ppt_thr = 0; % 0 mm = 0 inches
    if contains(current_var, "annual") % if the data is annual
        T_main{T_main_county_yr_idx,current_var} = sum(PRISM_yr_county_tdata(:,PRISM_key_idx) <= ppt_thr); % total annual data for currentcounty currentyear
    else % if the data is monthly
        for current_month = 1:length(month_names)
            if contains(current_var, month_names(current_month))
                break % all we need here is the current_month_idx
            end
        end
        month_idx = find(PRISM_month_idx == current_month);
        T_main{T_main_county_yr_idx,current_var} = sum(PRISM_yr_county_tdata(month_idx,PRISM_key_idx) <= ppt_thr); % total monthly data for currentcounty currentyear
    end
end
end

function [T_main] = T_main_assign_GDD(T_main,tabletype,month_names,T_main_county_yr_idx,PRISM_yr_county_tdata,PRISM_month_idx,i,j)
% insert special variables for no ppt days- both annual and monthly
nrows = length(T_main.Year);
if i == 1 && j == 1
    T_main_varname = horzcat('GDD_annual_',tabletype);
    T_main{:,T_main_varname} = zeros(nrows,1);
    for current_month = 1:length(month_names)
        T_main_varname = horzcat('GDD_',month_names{current_month},'_',tabletype);
        T_main{:,T_main_varname} = zeros(nrows,1);
    end
end

T_main_varnames = T_main.Properties.VariableNames;
var_startidx = find(contains(T_main_varnames, 'GDD_annual_')); % start index
var_endidx = length(T_main_varnames);

for var_idx = var_startidx:var_endidx
    current_var = T_main_varnames{var_idx};
    % Find GDD
    Tbase = 41; %%%%%%%%%%%%%%%%%% degrees farenheit
    PRISM_key = {'tmax','tmin'};
    PRISM_key_idx = [4,2]; % hard-coded:::::::: index of respective variables
    tmax = PRISM_yr_county_tdata(:,PRISM_key_idx(1));
    tmin = PRISM_yr_county_tdata(:,PRISM_key_idx(2));
    for day = 1:length(tmax)
        GDD(day) = (tmax(day) + tmin(day))/2 - Tbase;
    end

    if contains(current_var, "annual") % if the data is annual
        T_main{T_main_county_yr_idx,current_var} = sum(GDD>0); % total annual GDD for currentcounty currentyear
    else % if the data is monthly
        for current_month = 1:length(month_names)
            if contains(current_var, month_names(current_month))
                break % all we need here is the current_month_idx
            end
        end
        month_idx = find(PRISM_month_idx == current_month);
        T_main{T_main_county_yr_idx,current_var} = sum(GDD(month_idx)>0); % total monthly data for currentcounty currentyear
    end
end
end
