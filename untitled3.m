iMov=3;
for iP=1:3
subplot(3,1,iP);
cla;
plot(tempS(iChan,iMov).matsdf{iP}); hold on;
plot(S(iChan,iMov).mnsdf(:,iP), 'k-', 'LineWidth', 2);
title(sprintf('Movie %d', S(iChan,iMov).movID(iP)))
end

% Figures for methods
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

iChan=1;
iMov=1;

tAxis = [1:1:300000].*0.001; %0.001:0.001:300;
curRGR_resample = resample(S(iChan, iMov).mnsdf, 0.001*1000, 2.4*1000);
% convolve with gamma function
k = gampdf([-40:2.4:40],4,2);
curRGR = conv(curRGR_resample, k,'same');

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 300 800 500])
% subplot(3,1,1) % spike density function for each trial
% plot(tAxis, tempS(iChan,iMov).matsdf{1})
% xlabel('Time (s)')
% ylabel('Firing rate (spike/s)')
% box off
% title(sprintf('Cell %s, Movie %s',S(iChan, iMov).cellID, num2str(S(iChan, iMov).movID))) 

% subplot(3,1,2) % Averaged spike density function across trials (mean +- ste)
% plot(tAxis, S(iChan, iMov).mnsdf, 'k-', 'LineWidth', 2)
% hold on;
% tmpSte = std(tempS(iChan,iMov).matsdf{1}, 0, 2)./sqrt(size(tempS(iChan,iMov).matsdf{1},2));
% plot(tAxis, S(iChan, iMov).mnsdf+tmpSte, 'k-', 'LineWidth', 1)
% plot(tAxis, S(iChan, iMov).mnsdf-tmpSte, 'k-', 'LineWidth', 1)
% xlabel('Time (s)')
% ylabel('Firing rate (spike/s)')
% box off

subplot(3,1,3); cla; % MION-convolved regressor (resampled points are plotted with markers)
% plot(tAxis, S(iChan, iMov).rgrsMION, 'b-', 'LineWidth', 2)
hold on;
plot([1:125].*2.4, curRGR, 'bo', 'LineWidth', 2, 'MarkerFaceColor', 'w')
xlabel('Time (s)')
ylabel('arbitrary unit')
box off

print(gcf, '-depsc', fullfile(dirFig, ...
    sprintf('tor%smov%s_sdfRegressor', S(iChan, iMov).cellID, num2str(S(iChan, iMov).movID))));



% MION function
t = 0:2.4:300;
h = 16.4486 * ( -0.184/ 1.5 * exp(-t/ 1.5)...
+0.330/ 4.5 * exp(-t/ 4.5)...
+0.670/13.5 * exp(-t/13.5) );
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(h, 'LineWidth', 2)
box off
set(gca, 'XTick', [], 'YTick', [])
xlim([-3 125])
print(gcf, '-depsc', fullfile(dirFig, 'MION_IRF'));


% For each set of movies
for iMovSet = 1:3
    figure;
    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 300 400])
    
    subplot(3,1,1)
    plot(tAxis, tempS(iChan, (iMovSet-1)*3+1).matsdf{1});    
    hold on;
    plot(tAxis, S(iChan, (iMovSet-1)*3+1).mnsdf, 'k-', 'LineWidth', 2)
    title(sprintf('Movie %s', num2str(S(iChan, (iMovSet-1)*3+1).movID)))
    
    subplot(3,1,2)
    plot(tAxis, tempS(iChan, (iMovSet-1)*3+2).matsdf{1});
    hold on;
    plot(tAxis, S(iChan, (iMovSet-1)*3+2).mnsdf, 'k-', 'LineWidth', 2)
    title(sprintf('Movie %s', num2str(S(iChan, (iMovSet-1)*3+2).movID)))
    
    subplot(3,1,3)
    plot(tAxis, tempS(iChan, (iMovSet-1)*3+3).matsdf{1});
    hold on;
    plot(tAxis, S(iChan, (iMovSet-1)*3+3).mnsdf, 'k-', 'LineWidth', 2)
    title(sprintf('Movie %s', num2str(S(iChan, (iMovSet-1)*3+3).movID)))
    
    print(gcf, '-depsc', fullfile(dirFig, ...
    sprintf('tor%s_movieSet%d_SDF', S(iChan, iMov).cellID, iMovSet)));
end
    
    
iMov=3;
for iP=1:3
subplot(3,1,iP);
cla;
plot(tempS(iChan,iMov).matsdf{iP}); hold on;
plot(S(iChan,iMov).mnsdf(:,iP), 'k-', 'LineWidth', 2);
title(sprintf('Movie %d', S(iChan,iMov).movID(iP)))
end



figure;
set(gcf, 'Color', 'w', 'Position', [200 100 600 800], 'PaperPositionMode', 'auto')




figure
plot(S(iChan,iMov).rgrsMION_resample)
clf
subplot(2,1,1);
plot(S(iChan,iMov).rgrsMION_resample)
title(sprintf('Channel: %s, Movie: %d', S(iChan, iMov).cellID, S(iChan, iMov).movID(iMov)))
subplot(2,1,2);
dataBOLD
title(sprintf('Channel: %s, Movie: %s', S(iChan, iMov).cellID, num2str(S(iChan, iMov).movID)))

subplot(2,1,2)
title([])
title('fMRI time courses')

size(tempBOLD)

plot(tempBOLD(1,:))




plot(tempBOLD(iV,:), 'o-')
xlim([1 375])
help fft
help wavelet
help plotyy



tempBOLD = reshape(dataBOLD.catmvoltc{iMov}, nVox, nt);
[sortR, a] = sort(S(iChan, iMov).mapR(:), 'descend');
iV_pos=a(1);
iV_neg=a(end);
[~, iV_low] = min(abs(S(iChan,iMov).mapR(:))); % minimum (almost zero)

figure;
subplot(3,1,1); % positive correlation
ax=plotyy(1:375, S(iChan,iMov).rgrsMION_resample, 1:375, tempBOLD(iV_pos,:));
legend('Neural Regressor', 'fMRI Response')
title(sprintf('Channel: %s, Movie: %s, Voxel with positive correlation (R: %3.3f)',...
    S(iChan, iMov).cellID, num2str(S(iChan, iMov).movID), S(iChan,iMov).mapR(iV_pos))) 
axis(ax, 'tight')

subplot(3,1,2); % zero correlation
ax2=plotyy(1:375, S(iChan,iMov).rgrsMION_resample, 1:375, tempBOLD(iV_low,:));
title(sprintf('Channel: %s, Movie: %s, Voxel with no correlation (R: %3.3f)',...
    S(iChan, iMov).cellID, num2str(S(iChan, iMov).movID), S(iChan,iMov).mapR(iV_low)))
axis(ax2, 'tight')

subplot(3,1,3); % negative correlation
ax3=plotyy(1:375, S(iChan,iMov).rgrsMION_resample, 1:375, tempBOLD(iV_neg,:));
title(sprintf('Channel: %s, Movie: %s, Voxel with negative correlation (R: %3.3f)',...
    S(iChan, iMov).cellID, num2str(S(iChan, iMov).movID), S(iChan,iMov).mapR(iV_neg)))
axis(ax3, 'tight')


% Regressor per each movie
figCellResp = figure(101);
set(figCellResp, 'Color', 'w', 'Position', [200 100 540 770], 'PaperPositionMode', 'auto')
for iChan=1:length(S)
    figure(figCellResp)
    clf;
    
    for iMov=1:9
        figure(figCellResp);
        subplot(9,1,iMov);
        if isempty(S(iChan,iMov).rgrsMION_resample)
            title(sprintf('Movie %d', iMov))
            continue;
        end
        
        curCellID = S(iChan, iMov).cellID;
        tmpResample = resample(tempS(iChan,iMov).matsdf{1}, 0.001*1000, 2.4*1000);        
        tmpMn = mean(tmpResample,2);
        tmpSte = std(tmpResample, 0, 2)./sqrt(size(tmpResample,2));
%         errorbar([1:125].*2.4, tmpMn, tmpSte, 'LineWidth', 1.5)   
        plot([1:125].*2.4, tmpMn, 'LineWidth', 1.5)
        line([1:125; 1:125].*2.4, [tmpMn-tmpSte, tmpMn+tmpSte]', 'Color', 'b', 'LineWidth', 1.5)
%         plot([1:125].*2.4, S(iChan,iMov).rgrsMION_resample, 'LineWidth', 2)
        title(sprintf('Movie %d', S(iChan, iMov).movID))
        xlim([0 300])
        ylim([0 ceil(max(tmpMn+tmpSte))])
        set(gca, 'XMinorTick', 'on')
        
        if iMov==5
            ylabel('Firing rate (spikes/s)')
        end
        if iMov==9
            xlabel('Time (s)')
        end
    end
    
    fprintf('Cell %s (%d/%d)\n', curCellID, iChan, length(S))
    
    dirFig = '/einstein0/USRlab/projects/parksh/_labNote/_figs/';
    figName = sprintf('tor%s_avgSDF', curCellID);
    print(figCellResp, '-depsc', fullfile(dirFig,figName))
    
end





