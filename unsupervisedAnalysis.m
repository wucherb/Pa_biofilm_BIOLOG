% load data
tblCvOd = readtable("tidy_tables/tblCvOd.csv");

%% plot the distributions
figure(1)
subplot(2, 2, 1)
s = scatterhistogram(tblCvOd,'ODAero','CVAero', ...
    'HistogramDisplayStyle','smooth', ...
    'LineStyle','-');

subplot(2, 2, 2)
s = scatterhistogram(tblCvOd,'CVAnaero','CVAero', ...
    'HistogramDisplayStyle','smooth', ...
    'LineStyle','-');

subplot(2, 2, 3)
s = scatterhistogram(tblCvOd,'ODAnaero','CVAnaero', ...
    'HistogramDisplayStyle','smooth', ...
    'LineStyle','-');

subplot(2, 2, 4)
s = scatterhistogram(tblCvOd,'ODAnaero','ODAero', ...
    'HistogramDisplayStyle','smooth', ...
    'LineStyle','-');


%% plot the clustergram
% use this for only the nutrients with some evidence of growth
idx = (tblCvOd.odAnaeroPValue<0.05 | tblCvOd.odAeroPValue<0.05)& tblCvOd.log2fcCvAnaeroPValue < 0.05;

% use this to plot all 190 nutrients
%idx = 1:height(tblCvOd);

clustergram(tblCvOd{idx, {'ODAero','CVAero', 'ODAnaero','CVAnaero'}},  'ColumnPDist', 'correlation',...
    'ColumnLabels', {'ODAero','CVAero', 'ODAnaero','CVAnaero'},...
    'RowPDist', 'correlation', 'RowLabels', tblCvOd.Chemical(idx))

%% draw a pcoa that equivalent to the clustergram (using correlaiton distance)

% use this for only the nutrients with some evidence of growth
%idx = (tblCvOd.odAnaeroPValue<0.05 | tblCvOd.odAeroPValue<0.05)& tblCvOd.log2fcCvAnaeroPValue < 0.05;

% use this to plot all 190 nutrients
idx = 1:height(tblCvOd);

D = pdist(tblCvOd{idx, {'ODAero','CVAero', 'ODAnaero','CVAnaero'}}, 'correlation');
[cmdsScores, eigvals] = cmdscale(D, 2);


figure(2);
scatter(cmdsScores(:,1),cmdsScores(:,2), [],...
    tblCvOd.log2fcCvAnaero(idx),  'o', 'filled', 'MarkerEdgeColor','k')
xlabel('PCoA 1')
ylabel('PCoA 2')
grid on
axis equal square
title('PCoA based on correlation distance')
a = colorbar;
a.Label.String = 'log_2(FC) of biofilm in anerobic vs aerobic';
