% genFig_LN2018.m

% compare Avalanche 069a Matcha 097a time courses (these two show really
% similar maps)
% they are in _orgData


iMov = 1;
matCell = '081b'; %'097a';
avaCell = '117a'; %'069a';
load(sprintf('/procdata/parksh/Mat/mmov%dsig%s.mat', iMov, matCell));
dat_mat = dat.s;
load(sprintf('/procdata/parksh/Ava/amov%dsig%s.mat', iMov, avaCell));
dat_ava = dat.s;

figure;
subplot(2,1,1)
for iTrial = 1:length(dat_ava)
line(repmat(dat_ava{iTrial}', 2, 1), [ones(size(dat_ava{iTrial})).*iTrial ones(size(dat_ava{iTrial})).*(iTrial+1)]', 'Color', 'k')
hold on
end
title(sprintf('Avalanche cell%s', avaCell))
subplot(2,1,2)
for iTrial = 1:length(dat_mat)
line(repmat(dat_mat{iTrial}', 2, 1), [ones(size(dat_mat{iTrial})).*iTrial ones(size(dat_mat{iTrial})).*(iTrial+1)]', 'Color', 'k')
hold on
end
title(sprintf('Matcha cell%s', matCell))

iMov = 1;
S_1 = createCellRegressor_indMov('/procdata/parksh/Dan', {'18', '23', '28'}, iMov, 1000); % last input is sigma (for SDF computation) in msec
S_2 = createCellRegressor_indMov('/procdata/parksh/Ava/_beforeSelection', {'117a'}, iMov, 1000); %
% S_2 = createCellRegressor_indMov('/procdata/parksh/Spi/2018Jan_movie', {'55'}, iMov, 1000); %%'/procdata/parksh/Spi/2018Jan_movie
tempC = bone(15); % jet(size(S(1).matsdf{1}, 2))*0.75; % temporary ColorOrder for line plot

figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 600])
% case of two example neurons
SP(1) = subplot('position', [0.1 0.55 0.8 0.35]);
SP(2) = subplot('position', [0.1 0.15 0.8 0.35]);
set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(SP(1), S_1.matsdf{1}, 'LineWidth', 2)
plot(SP(2), S_2.matsdf{1}, 'LineWidth', 2)
set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
% set(SP, 'YLim', [0 30], 'YTick', [0 30], 'YTickLabel', [0 30], 'XTickLabel', [])

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 853 210])
plot(zscore(S_2.mnsdf), 'm-')
hold on
plot(zscore(S_1.mnsdf))


nameSubj = 'Dav';
iMov = 4;
setCellIDs = {'03', '17', '19'};
nCell=length(setCellIDs);

S = createCellRegressor_indMov(fullfile('/procdata/parksh', nameSubj), setCellIDs, iMov, 1000); 

figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 800 860])
% case of five example neurons
startY = 0.07;
height = 0.15;
gap = 0.03;
for iCell = 1:nCell
    SP(iCell) = subplot('position', [0.1 startY+height*(iCell-1)+gap*(iCell-1) 0.8 height]);
end
set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
for iCell = 1:nCell
    plot(SP(iCell), S(iCell).matsdf{1}, 'LineWidth', 2)
end

set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
% set(SP(3:5), 'YLim', [0 30], 'YTick', [0 30], 'YTickLabel', [0 30])
% set(SP(1:2), 'YLim', [0 10], 'YTick', [0 10], 'YTickLabel', [0 10])
set(SP, 'XTickLabel', [])
% set(SP(1), 'XTickLabel', 0:50:300)
xlabel(SP(1), 'Time (sec)')
ylabel(SP(3), 'Spike / sec')
% save
print(figSDF, fullfile(dirFig, sprintf('exampleSDF_%s_cell%s_trialbytrial', nameSubj, cell2mat(setCellIDs))), '-depsc')

% FR_dT_mat = createCellRegressor_indMov_discreteTime('/procdata/parksh/Mat', {'097a'}, 4, 0.1); % 
% FR_dT_ava = createCellRegressor_indMov_discreteTime('/procdata/parksh/Ava/_beforeSelection', {'069a'}, 4, 0.1); % 

iMov = 6;
FR_dT_1 = createCellRegressor_indMov_discreteTime('/procdata/parksh/Mat', {'081b'}, iMov, 0.1);
FR_dT_2 = createCellRegressor_indMov_discreteTime('/procdata/parksh/Ava/_beforeSelection', {'117a'}, iMov, 0.1); %
% FR_dT_2 = createCellRegressor_indMov_discreteTime('/procdata/parksh/Spi/2018Jan_movie', {'55'}, iMov, 0.1); %
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 853 327])
subplot(2,1,1)
imagesc(FR_dT_1.matFR{1}')
set(gca, 'CLim', [0 1])
colormap(flipud(gray))
subplot(2,1,2)
imagesc(FR_dT_2.matFR{1}')
set(gca, 'CLim', [0 1])
