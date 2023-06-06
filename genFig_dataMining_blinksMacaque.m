load /nifvault/procdata/parksh/_macaque/_Harish/matSDF_forshp.mat;
load /nifvault/procdata/parksh/_macaque/_Harish/ed_forshp.mat;

set(0, 'defaultfigurecolor', [1 1 1]);
% trying most conservative ways to see the difference between pre - during
% - post blink period in cell activity (firing rate)

% skip subject #7 because the pupil data range is different and detection
% of blinks based on thresholding does not work 
iSubj = 9;
iMov = 3;
iTrial = 1;
iCell = 1;

figure;
subplot(211);
plot(ed(iSubj).mov(iMov).isblink(:,iTrial),'k');
hold on;
plot(ed(iSubj).mov(iMov).rawpd(:,iTrial),'.-g');
% plot(ed(iSubj).mov(iMov).pd(:,iTrial),'r');
legend('blink', 'rawPD') %, 'PD')
subplot(212);
triald = matSDF(iSubj).FR_30fps(iCell,iMov).matFR{iTrial};
plot(triald(1:1000,1));
legend('trial FR');
title(sprintf('Monkey %d Movie %d Trial %d Cell %d', iSubj, iMov, iTrial, iCell));

% identify the blinks more conservative way

% sub-select the blinks containing pupil signal below -10

for iSubj = 7:9
    for iMov = 1:3
        cMap = pink(size(ed(iSubj).mov(iMov).rawpd, 2));
        figure(2)
        clf
        for iTrial = 1:size(ed(iSubj).mov(iMov).rawpd, 2)
            histogram(ed(iSubj).mov(iMov).rawpd(:,iTrial), -12:0.5:12, 'FaceColor', cMap(iTrial,:));
            hold on;
        end
        title(sprintf('Distribution of analog pupil signal: monkey %d movie %d', iSubj, iMov));
        
        input('')
    end
end



%%
iSubj = 9; %7; %:9
iMov = 1; %:3
for iTrial = [1 2 5] %1:size(ed(iSubj).mov(iMov).rawpd, 2)
    
    temp = diff(ed(iSubj).mov(iMov).isblink(:,iTrial)); % 1 for blink start, we will loop for 1s in this vector
    curLocBlinkOnset = find(temp>0)+1;
    locBlinkOffset = find(temp<0)+1;
    if locBlinkOffset(1)<curLocBlinkOnset(1)
        locBlinkOffset = locBlinkOffset(2:end);
    elseif curLocBlinkOnset(end)>locBlinkOffset(end)
        curLocBlinkOnset = curLocBlinkOnset(1:end-1);
    end
    durBlink = [];
    durBlink_frame = diff([curLocBlinkOnset locBlinkOffset], 1, 2);
    durBlink_ms = diff([curLocBlinkOnset locBlinkOffset], 1, 2).*33;
         
    minPD_blink=[]; tPD_blink={};
    for iBlink = 1:length(curLocBlinkOnset) % we want to compute the mean of the PD during blinks to subselect blinks
        tPD_blink{iBlink} = ed(iSubj).mov(iMov).rawpd(curLocBlinkOnset(iBlink):locBlinkOffset(iBlink), iTrial);
        minPD_blink(iBlink,1) = min(tPD_blink{iBlink});
    
%         figure(3);cla;
%         plot(tempPD)
%         input('')
    end
    
    indValidBlink=[];
    indValidBlink = find(minPD_blink < -5);
    if isempty(indValidBlink)
        continue;
    end
    
    resultsBlink(iTrial).durBlink_ms = durBlink_ms(indValidBlink);
    resultsBlink(iTrial).durBlink_frame = durBlink_frame(indValidBlink);
    resultsBlink(iTrial).locBlinkOnset = curLocBlinkOnset(indValidBlink);
    resultsBlink(iTrial).locBlinkOffset = locBlinkOffset(indValidBlink);
    resultsBlink(iTrial).tPD_blink = tPD_blink(indValidBlink);
    
%     [a, indSort] = sort(minPD_blink, 'ascend'); % from smallest voltage to largest 
%     figure(100); clf;
%     for i = 1:length(indValidBlink)
%         figure(100)
%         plot(tPD_blink{indValidBlink(i)});
%         hold on;
%     end
%     title(sprintf('Monkey %d Movie %d Trial %d', iSubj, iMov, iTrial));
%     input('')
    
%     rawpd_blink{iTrial} = ed(iSubj).mov(iMov).rawpd(find(ed(iSubj).mov(iMov).isblink(:,iTrial)>0), iTrial);
%     
end

iCell = 1;
fig_cellTS_blink = figure;
% fig_cellTS_blinkpre = figure;
% fig_cellTS_blinkpost = figure;

setTrial = [1 2 5];
for iT = 1:length(setTrial)
    
    iTrial = setTrial(iT);
    
    curLocBlinkOnset = resultsBlink(iTrial).locBlinkOnset;
    curLocBlinkOffset = resultsBlink(iTrial).locBlinkOffset;
    curDurBlink_frame = resultsBlink(iTrial).durBlink_frame;
        
    figure(fig_cellTS_blink);
    subplot(1,3,iT);
    for iBlink = 1:length(curLocBlinkOnset)
        
        tS_blink{iBlink} = matSDF(iSubj).FR_30fps(iCell, iMov).matFR{1}(curLocBlinkOnset(iBlink):curLocBlinkOffset(iBlink), iTrial);
        tS_blink_pre{iBlink} = matSDF(iSubj).FR_30fps(iCell, iMov).matFR{1}(curLocBlinkOnset(iBlink)-curDurBlink_frame(iBlink)-1:curLocBlinkOnset(iBlink)-1, iTrial);
        tS_blink_post{iBlink} = matSDF(iSubj).FR_30fps(iCell, iMov).matFR{1}(curLocBlinkOffset(iBlink)+1:curLocBlinkOffset(iBlink)+1+curDurBlink_frame(iBlink), iTrial);
        
        plot(tS_blink_pre{iBlink}); hold on;
    end
    title(sprintf('Monkey %d Movie %d Trial %d Cell %d', iSubj, iMov, iTrial, iCell));

    tS_mean_blink = cellfun(@mean, tS_blink);
    tS_mean_blink_pre = cellfun(@mean, tS_blink_pre);
    tS_mean_blink_post = cellfun(@mean, tS_blink_post);
    
    resultsBlink(iTrial).tS_blink = tS_blink;
    resultsBlink(iTrial).tS_blink_pre = tS_blink_pre;
    resultsBlink(iTrial).tS_blink_post = tS_blink_post;
    resultsBlink(iTrial).tS_mean_blink = tS_mean_blink;
    resultsBlink(iTrial).tS_mean_blink_pre = tS_mean_blink_pre;
    resultsBlink(iTrial).tS_mean_blink_post = tS_mean_blink_post;
end
    
catMean_pre = cat(2, resultsBlink(:).tS_mean_blink_pre)'; 
catMean_post = cat(2, resultsBlink(:).tS_mean_blink_post)'; 
catMean_during = cat(2, resultsBlink(:).tS_mean_blink)'; 

matMeanTS = cat(2, catMean_pre, catMean_post, catMean_during);
R = corr(matMeanTS, 'type', 'Spearman');

figure
subplot(1,3,1)
plot(catMean_pre, catMean_during, 'o')
xlabel('before')
ylabel('during')
subplot(1,3,2)
plot(catMean_post, catMean_during, 'o')
xlabel('after')
ylabel('during')
subplot(1,3,3)
plot(catMean_pre, catMean_post, 'o')
xlabel('before')
ylabel('after')


figure;
plot(ed(iSubj).mov(iMov).isblink(:,iTrial),'k');
hold on;
plot(ed(iSubj).mov(iMov).rawpd(:,iTrial),'.-g');
% plot(ed(iSubj).mov(iMov).pd(:,iTrial),'r');
legend('blink', 'rawPD') %, 'PD')

title(sprintf('Monkey %d Movie %d Trial %d Cell %d', iSubj, iMov, iTrial, iCell));