% seeCorrMap_tSeries.m

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
% BlockANa corrmaps should be already loaded in the workplace
global DATA DSP

z_timecourse = (DSP.timecourse-nanmean(DSP.timecourse))./nanstd(DSP.timecourse);
z_currgr     = (DATA.stim.curstimmodel-nanmean(DATA.stim.curstimmodel))./nanstd(DATA.stim.curstimmodel);

figure;
set(gcf, 'Color', 'k', 'PaperPositionMode', 'auto', 'Position', [100 100 970 230])
plot(DATA.ftimes, z_currgr, 'm-', 'LineWidth', 2)

hold on
plot(DATA.ftimes, z_timecourse, 'w-', 'LineWidth', 2)
hold on

set(gca, 'color', 'k', 'TickDir', 'out', 'box', 'off')
set(gca, 'LineWidth', 2)
set(gca, 'XColor', 'w')
set(gca, 'YColor', 'w')
set(gca, 'FontSize', 15)
set(gca, 'FontSize', 12)

figName = sprintf('exampleCorrTseries_highGammaBLP_vox_%s', strrep(num2str(DSP.or.activevoxel), ' ', '_'));
print(gcf, fullfile(dirFig, figName), '-r300', '-dtiff')
