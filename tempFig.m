figTS = figure;
set(figTS, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 295])
clf;
plot(DATA.ftimes, z_timecourse, 'b-', 'LineWidth', 2)
hold on
plot(DATA.ftimes, z_currgr, 'm-', 'LineWidth', 2)
set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 15)
voxi = DSP.or.activevoxel;
arr = DSP.proc.scalarmap_3d;
title(sprintf('Vox Coords = [%d %d %d], r = %.2f', voxi(1), voxi(2), voxi(3), arr(voxi(1),voxi(2),voxi(3))))


figXCorr = figure;
set(figXCorr, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 400 300])
clf;
timecourse = squeeze(curFMRITC(voxi(1),voxi(2),voxi(3),:));
mntc = nanmean(timecourse);
timecourse(isnan(timecourse)) = mntc;

% -1 below for mion
[xc,lags] = xcorr(zscore(curRGR),-1*zscore(timecourse),20,'coeff');
[TR,NR] = getTRandNR();
lags = lags*TR;
plot(lags, xc)

