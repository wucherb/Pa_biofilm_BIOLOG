% load the xlsx file
tblXls = readtable('biolog_data/pabiolog_master_update.xlsx');


%% create Nutrient/Class table
tblNutrientClass = tblXls(:, {'Chemical' 'SeriesMoA'});
tblNutrientClass = unique(tblNutrientClass, 'rows');
% save the tidy table
writetable(tblNutrientClass, 'tidy_tables/tblNutrientClass.csv', 'Delimiter',',');

%% display a ranked list of nutrient classes
tblResponses = tblXls(:, {'Chemical', 'replicate', 'Plate', 'CV_a', 'Odb_a', 'CV_an', 'Odb_an'});
tblResponses = stack(tblResponses, 4:7, 'NewDataVariableName','value', 'IndexVariableName','type_aero');
% split the type_aero into type and aero
for i = 1:height(tblResponses)
    straux = strsplit(string(tblResponses.type_aero(i)), '_');
    tblResponses.type{i} = straux{1};
    tblResponses.aero{i} = straux{2};
end
% convert aero to boolean (0 for anaerobic, 1 for aerobic)
tblResponses.aero = strcmp(tblResponses.aero, 'a');
tblResponses.type_aero = [];

%% Compute the CV and OD values corrected for negative control
idxAero = tblResponses.aero;

% AEROBIC
tblCV = tblResponses(strcmp(tblResponses.type, 'CV') & idxAero, :);
tblOD = tblResponses(strcmp(tblResponses.type, 'Odb') & idxAero, :);

% OD aerobic
mdlOD = fitlme(tblOD, 'value ~ Chemical + (1|Plate:replicate)');
tblOdAero = dataset2table(mdlOD.Coefficients);
tblOdAero = tblOdAero(:, {'Name' 'Estimate' 'pValue'});
tblOdAero.Properties.VariableNames = {'Chemical' 'ODAero' 'odAeroPValue'};
% CV aerobic
mdlCv = fitlme(tblCV, 'value ~ Chemical + (1|Plate:replicate)');
tblCvAero = dataset2table(mdlCv.Coefficients);
tblCvAero = tblCvAero(:, {'Name' 'Estimate' 'pValue'});
tblCvAero.Properties.VariableNames = {'Chemical' 'CVAero' 'cvAeroPValue'};
% remove the row corresponding to the intercept
tblCvOdAero = innerjoin(tblCvAero, tblOdAero);
tblCvOdAero(1, :) = [];

% ANEROBIC
tblCV = tblResponses(strcmp(tblResponses.type, 'CV') & ~idxAero, :);
tblOD = tblResponses(strcmp(tblResponses.type, 'Odb') & ~idxAero, :);

% OD anaerobic
mdlOD = fitlme(tblOD, 'value ~ Chemical + (1|Plate:replicate)');
tblOdAnaero = dataset2table(mdlOD.Coefficients);
tblOdAnaero = tblOdAnaero(:, {'Name' 'Estimate' 'pValue'});
tblOdAnaero.Properties.VariableNames = {'Chemical' 'ODAnaero' 'odAnaeroPValue'};
% CV anaerobic
mdlCv = fitlme(tblCV, 'value ~ Chemical + (1|Plate:replicate)');
tblCvAnaero = dataset2table(mdlCv.Coefficients);
tblCvAnaero = tblCvAnaero(:, {'Name' 'Estimate' 'pValue'});
tblCvAnaero.Properties.VariableNames = {'Chemical' 'CVAnaero' 'cvAnaeroPValue'};
% remove the row corresponding to the intercept
tblCvOdAnaero = innerjoin(tblCvAnaero, tblOdAnaero);
tblCvOdAnaero(1, :) = [];

% compile aerobic and anaerobic data
tblCvOd = innerjoin(tblCvOdAero, tblCvOdAnaero);

%% remove 'Chemical_' from the chemical name
tblCvOd.Chemical = strrep(tblCvOd.Chemical, 'Chemical_', '');

%% fcreate auxiliary variables to determine fold changes
tblResponses.logValue = log(tblResponses.value);
tblResponses.type = string(tblResponses.type);
tblResponses.anaero = ~(tblResponses.aero);

%% fitlm to determine the change in biofilm formation when switching to anaerobic
mdlCvFc = fitlme(tblResponses, 'logValue ~ anaero*Chemical + (1|Plate:replicate)',...
    'Exclude',~strcmp(tblResponses.type, 'CV'));

% Keep only the fold changes in anaerobic
tblCvFoldChange = dataset2table(mdlCvFc.Coefficients);
% convert from natural base to log2FC
tblCvFoldChange.Estimate = tblCvFoldChange.Estimate / log(2);
% clean up column names
tblCvFoldChange.Properties.VariableNames{1} = 'Chemical';
idx = contains(tblCvFoldChange.Chemical, ':anaero_1');
tblCvFoldChange(~idx, :) = [];
tblCvFoldChange.Chemical = strrep(tblCvFoldChange.Chemical, 'Chemical_', '');
tblCvFoldChange.Chemical = strrep(tblCvFoldChange.Chemical, ':anaero_1', '');

% add these values to the tblCvOd 
tblCvFoldChange = tblCvFoldChange(:, {'Chemical' 'Estimate' 'pValue'});
tblCvFoldChange.Properties.VariableNames = {'Chemical' 'log2fcCvAnaero' 'log2fcCvAnaeroPValue'};
tblCvOd = innerjoin(tblCvOd, tblCvFoldChange);


%% fitlm to determine the change in OD when switching to anaerobic
mdlOdFc = fitlme(tblResponses, 'logValue ~ anaero*Chemical + (1|Plate:replicate)',...
    'Exclude', strcmp(tblResponses.type, 'CV'));

% Keep only the fold changes in anaerobic
tblOdFoldChange = dataset2table(mdlOdFc.Coefficients);
% convert from natural base to log2FC
tblOdFoldChange.Estimate = tblOdFoldChange.Estimate / log(2);
% clean up column names
tblOdFoldChange.Properties.VariableNames{1} = 'Chemical';
idx = contains(tblOdFoldChange.Chemical, ':anaero_1');
tblOdFoldChange(~idx, :) = [];
tblOdFoldChange.Chemical = strrep(tblOdFoldChange.Chemical, 'Chemical_', '');
tblOdFoldChange.Chemical = strrep(tblOdFoldChange.Chemical, ':anaero_1', '');

% add these values to the tblCvOd 
tblOdFoldChange = tblOdFoldChange(:, {'Chemical' 'Estimate' 'pValue'});
tblOdFoldChange.Properties.VariableNames = {'Chemical' 'log2fcOdAnaero' 'log2fcOdAnaeroPValue'};
tblCvOd = innerjoin(tblCvOd, tblOdFoldChange);

%% save the table
writetable(tblCvOd, 'tidy_tables/tblCvOd.csv', 'Delimiter',',');

%%
tblCvOd = sortrows(tblCvOd, 'log2fcCvAnaero')
idx = (tblCvOd.odAnaeroPValue<0.05 | tblCvOd.odAeroPValue<0.05)& tblCvOd.log2fcCvAnaeroPValue < 0.05;

figure(5)
barh(tblCvOd.log2fcCvAnaero(idx)/log(2))
grid on
set(gca, 'YTick', 1:sum(idx), 'YTickLabel', tblCvOd.Chemical(idx), 'TickLabelInterpreter', 'none')
xlabel('log_2(FC) of biofilm in anaerobic vs aerobic');
title({'Impact of anaerobiosis on biofilm formation'...
    'for nutrients with evidence they can be utilized'});

