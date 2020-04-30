function [] = seeCorrMap(S, dataBOLD, iChan, iMov)

global STDPATH DSP DATA GH tempS

% change the frac var map
DSP.proc.fracvarmap_3d = S(iChan, iMov).mapR;

% Trial-by-trial SDF
curmatsdf = tempS(iChan,iMov).matsdf{mod(iMov-1, 3)+1};

figure(100);

subplot('Position',[0.05 0.1 0.45 0.4]);
plot(curmatsdf)
ylabel('Firing rate (spikes/s)')
xlabel('Time (ms)')

% Averaged SDF across trials
subplot('Position',[0.05 0.53 0.45 0.2]);
plot(S(iChan, iMov).mnsdf, 'k-', 'LineWidth', 2)

% MION convolved
subplot('Position',[0.05 0.76 0.45 0.2]);
plot(S(iChan, iMov).rgrsMION, 'b-', 'LineWidth', 2)

% show time course of selected voxel
DSP.proc.fmri_tc_3d= dataBOLD.catsa5vmvoltc{iMov};

subplot('Position',[0.53 0.53 0.45 0.43]);
axx=plotyy(DATA.ftimes, S(iChan, iMov).rgrsMION_resample, DATA.ftimes, DSP.timecourse);
axis(axx, 'tight')
legend('Neural regressor', 'Voxel time course')
xlabel('Time (s)')


% plot(DATA.ftimes, DSP.timecourse,'b-','LineWidth',2)




