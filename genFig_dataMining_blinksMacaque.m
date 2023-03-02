load matSDF_forshp.mat;
load ed_forshp.mat;

iSubj = 8;
iMov = 1;
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