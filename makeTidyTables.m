% load the xlsx file
tblXls = readtable('biolog_data/pabiolog_master.xlsx');


%% plot the joint distribution of CVs and ODs
tblXls.negativeControl = strcmp(tblXls.Chemical, 'NegativeControl');
figure(1)
subplot(2, 1, 1)
s = scatterhistogram(tblXls,'Odb_an','CV_an', ...
    'GroupVariable','negativeControl', 'HistogramDisplayStyle','smooth', ...
    'LineStyle','-', 'XLimits',[0.05 1.1], 'YLimits',[0.35 2.4]);
title('Anaerobic')

subplot(2, 1, 2)
s = scatterhistogram(tblXls,'Odb_a','CV_a', ...
    'GroupVariable','negativeControl', 'HistogramDisplayStyle','smooth', ...
    'LineStyle','-', 'XLimits',[0.05 1.1], 'YLimits',[0.35 2.4]);
title('Aerobic')

%% create Nutrient/Class table
tblNutrientClass = tblXls(:, {'Chemical' 'SeriesMoA'});
tblNutrientClass = unique(tblNutrientClass, 'rows');
% save the tidy table
writetable(tblNutrientClass, 'tidy_tables/tblNutrientClass.csv', 'Delimiter',',');

%% display a ranked list of nutrient classes
sortrows(groupcounts(tblNutrientClass, 'SeriesMoA'), 'GroupCount', 'descend')


tblResponses = tblXls(:, {'Chemical', 'CV_a', 'Odb_a', 'CV_an', 'Odb_an'});
tblResponses = stack(tblResponses, 2:5, 'NewDataVariableName','value', 'IndexVariableName','type_aero');
% split the type_aero into type and aero
for i = 1:height(tblResponses)
    straux = strsplit(string(tblResponses.type_aero(i)), '_');
    tblResponses.type{i} = straux{1};
    tblResponses.aero{i} = straux{2};
end
% convert aero to boolean (0 for anaerobic, 1 for aerobic)
tblResponses.aero = strcmp(tblResponses.aero, 'a');
tblResponses.type_aero = [];

%% plot the distribution of CV values
idx = strcmp(tblResponses.Chemical, 'NegativeControl');
figure(2)
subplot(2, 1, 1)
histogram(tblResponses.value(strcmp(tblResponses.type, 'CV') & ~idx, :));
hold on
histogram(tblResponses.value(strcmp(tblResponses.type, 'CV') & idx, :));
hold off
xlabel('CV')
ylabel('all data')

subplot(2, 1, 2)
histogram(tblResponses.value(strcmp(tblResponses.type, 'Odb') & ~idx, :));
hold on
histogram(tblResponses.value(strcmp(tblResponses.type, 'Odb') & idx, :));
hold off
xlabel('OD')
ylabel('all data')

%% Compute the CV and OD values corrected for negative control
idxAero = tblResponses.aero;
tblCV = tblResponses(strcmp(tblResponses.type, 'CV') & idxAero, :);
tblOD = tblResponses(strcmp(tblResponses.type, 'Odb') & idxAero, :);

mdlCv = fitlm(tblCV, 'value ~ Chemical');
mdlOD = fitlm(tblOD, 'value ~ Chemical');


tblCvAero = mdlCv.Coefficients;
tblCvAero.SE = []; tblCvAero.tStat = [];
tblCvAero.Properties.VariableNames = {'CVAero' 'cvAeroPValue'};
tblCvAero.Chemical = tblCvAero.Properties.RowNames;

tblOdAero = mdlOD.Coefficients;
tblOdAero.SE = []; tblOdAero.tStat = [];
tblOdAero.Properties.VariableNames = {'ODAero' 'odAeroPValue'};
tblOdAero.Chemical = tblOdAero.Properties.RowNames;

tblCvOdAero = innerjoin(tblCvAero, tblOdAero);
tblCvOdAero(1, :) = [];


tblCV = tblResponses(strcmp(tblResponses.type, 'CV') & ~idxAero, :);
tblOD = tblResponses(strcmp(tblResponses.type, 'Odb') & ~idxAero, :);

mdlCv = fitlm(tblCV, 'value ~ Chemical');
mdlOD = fitlm(tblOD, 'value ~ Chemical');


tblCvAnaero = mdlCv.Coefficients;
tblCvAnaero.SE = []; tblCvAnaero.tStat = [];
tblCvAnaero.Properties.VariableNames = {'CVAnaero' 'cvAnaeroPValue'};
tblCvAnaero.Chemical = tblCvAnaero.Properties.RowNames;

tblOdAnaero = mdlOD.Coefficients;
tblOdAnaero.SE = []; tblOdAnaero.tStat = [];
tblOdAnaero.Properties.VariableNames = {'ODAnaero' 'odAnaeroPValue'};
tblOdAnaero.Chemical = tblOdAnaero.Properties.RowNames;

tblCvOdAnaero = innerjoin(tblCvAnaero, tblOdAnaero);
tblCvOdAnaero(1, :) = [];
tblCvOd = innerjoin(tblCvOdAero, tblCvOdAnaero);

%% change the order of columns so Chemical comes first
tblCvOd = tblCvOd(:, [3 1:2 4:end])
tblCvOd.Chemical = strrep(tblCvOd.Chemical, 'Chemical_', '');

%% plot the joint distribution of corrected averaged CVs and ODs

figure(3)
subplot(2, 1, 1)
s = scatterhistogram(tblCvOd,'ODAnaero','CVAnaero', ...
    'HistogramDisplayStyle','smooth', ...
    'LineStyle','-', 'XLimits',[-0.1 0.8], 'YLimits',[-0.4 1.5]);
title('Anaerobic')

subplot(2, 1, 2)
s = scatterhistogram(tblCvOd,'ODAero','CVAero', ...
    'HistogramDisplayStyle','smooth', ...
    'LineStyle','-', 'XLimits',[-0.1 0.8], 'YLimits',[-0.4 1.5]);
title('Aerobic')

%% color by pvalue
figure(4)
subplot(2, 2, 1)
gscatter(tblCvOd.ODAnaero, tblCvOd.CVAnaero, tblCvOd.odAnaeroPValue<0.05)
xlabel('ODAnaero')
ylabel('CVAnaero')
grid on
legend('n.s.', 'P<0.05')
title('P-value for OD')

subplot(2, 2, 2)
gscatter(tblCvOd.ODAnaero, tblCvOd.CVAnaero, tblCvOd.cvAnaeroPValue<0.05)
xlabel('ODAnaero')
ylabel('CVAnaero')
grid on
legend('n.s.', 'P<0.05')
title('P-value for CV')

subplot(2, 2, 3)
gscatter(tblCvOd.ODAero, tblCvOd.CVAero, tblCvOd.odAeroPValue<0.05)
xlabel('ODAero')
ylabel('CVAero')
grid on
legend('n.s.', 'P<0.05')
title('P-value for OD')

subplot(2, 2, 4)
gscatter(tblCvOd.ODAero, tblCvOd.CVAero, tblCvOd.cvAeroPValue<0.05)
xlabel('ODAero')
ylabel('CVAero')
grid on
legend('n.s.', 'P<0.05')

%% fitlm to determine the change in biofilm formation when switching to anaerobic
tblResponses.logValue = log(tblResponses.value);
tblResponses.type = string(tblResponses.type);
tblResponses.anaero = ~(tblResponses.aero);
%tblResponses.type = categorical(tblResponses.type);
%tblResponses.type = setcats(tblResponses.type, {'Odb' 'CV'})
mdlScaled = fitlm(tblResponses, 'logValue ~ anaero*Chemical',...
    'Exclude',~strcmp(tblResponses.type, 'CV'));

%% plot the nutrients that induce biggest change in CV
tblCvChange = mdlScaled.Coefficients;
% keep only the rows reporting CV changes in anaerobic
idx = contains(tblCvChange.Properties.RowNames, ':anaero_1');
tblCvChange(~idx, :) = [];
tblCvChange.Chemical = strrep(tblCvChange.Properties.RowNames, 'Chemical_', '');
tblCvChange.Chemical = strrep(tblCvChange.Chemical, ':anaero_1', '');

%% add these values to the tblCvOd 
tblCvChange = tblCvChange(:, {'Chemical' 'Estimate' 'pValue'});
tblCvChange.Properties.VariableNames = {'Chemical' 'log2fcCvAnaero' 'log2fcCvAnaeroPValue'};
tblCvOd = innerjoin(tblCvOd, tblCvChange)
tblCvOd = sortrows(tblCvOd, 'log2fcCvAnaero');

%% keep only the ones there is evidence they can utilize the nutrient
% bar plot of the change in biofilm formation when switching from aerobic
% to anaerobic conditions
% but only for nutrients where there is evidence they can be utilized for
% growth above engative control, either aerobically or anaerobically.
idx = (tblCvOd.odAnaeroPValue<0.05 | tblCvOd.odAeroPValue<0.05)& tblCvOd.log2fcCvAnaeroPValue < 0.05;

figure(5)
barh(tblCvOd.log2fcCvAnaero(idx)/log(2))
grid on
set(gca, 'YTick', 1:sum(idx), 'YTickLabel', tblCvOd.Chemical(idx), 'TickLabelInterpreter', 'none')
xlabel('log_2(FC) of biofilm in anaerobic vs aerobic');
title({'Impact of anaerobiosis on biofilm formation'...
    'for nutrients with evindence they can be utilized'});

%% save the table
writetable(tblCvOd, 'tidy_tables/tblCvOd.csv', 'Delimiter',',');


