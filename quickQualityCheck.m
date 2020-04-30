% quickQualityCheck.m
%
% Script for quick quality check of the data
% Just gathered common lines that I use frequently 
% to check the quality of data and get some hint of what's going on
% Usually consistenty acorss trials, plot each trial and averaged one, etc.


% trial-to-trial BLP time series for each band
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto'); clf;
for iBP = 1:4
    tempBLP = squeeze(dat.blp(:,1:200:60200, iBP));
subplot(4,1,iBP)
plot(tempBLP') % a little bit downsampling
axis tight
hold on;
title(sprintf('%d Hz - %d Hz', dat.h.blpfreq(iBP,1), dat.h.blpfreq(iBP,2)))

% plot(mean(tempBLP), 'k-', 'LineWidth', 3)
% errorbar(mean(tempBLP), (std(tempBLP).*ones(1, size(tempBLP,2)))./sqrt(size(tempBLP,2)-1),'k-')

end

% trial-to-trial BLP time series for each band (z-scored)
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto'); clf;
for iBP = 1:4
    matBLP = squeeze(dat.blp(:,:, iBP));
    zscoreBLP = zscore(matBLP, 0, 2);
    
    subplot(4,1,iBP)
    plot(zscoreBLP') % a little bit downsampling
    axis tight
    hold on;
    line([0 60200], [5 5], 'Color', 'k', 'LineStyle', '--') % line at 5*std
    
    title(sprintf('%d Hz - %d Hz', dat.h.blpfreq(iBP,1), dat.h.blpfreq(iBP,2)))
    
    % plot(mean(tempBLP), 'k-', 'LineWidth', 3)
    
    % errorbar(mean(tempBLP), (std(tempBLP).*ones(1, size(tempBLP,2)))./sqrt(size(tempBLP,2)-1),...
    %     'k-')

end

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto'); clf;
for iBP=1:4
subplot(4,1,iBP)
plot(BLPRGR(iMovie).matBLP{iBP})
hold on
plot(BLPRGR(iMovie).meanBLP{iBP}, 'k','LineWidth', 2)
title(sprintf('%d Hz - %d Hz', BLPRGR(1).blpfreq(iBP,1), BLPRGR(1).blpfreq(iBP,2)))
end

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto'); clf;
for iBP = 1:4
    tempBLP = squeeze(dat.blp(:,1:200:60200, iBP));
    meanBLP = mean(tempBLP);
    steBLP = (std(tempBLP).*ones(1, size(tempBLP,2)))./sqrt(size(tempBLP,2)-1);
subplot(4,1,iBP)
plot(meanBLP, 'k-') 
hold on;
line([1:size(tempBLP,2); 1:size(tempBLP,2)], [meanBLP+steBLP; meanBLP-steBLP], 'Color', 'k')
xlim([1 length(meanBLP)])
% axis tight

title(sprintf('%d Hz - %d Hz', dat.h.blpfreq(iBP,1), dat.h.blpfreq(iBP,2)))


% errorbar(mean(tempBLP), (std(tempBLP).*ones(1, size(tempBLP,2)))./sqrt(size(tempBLP,2)-1),...
%     'k-')

end

