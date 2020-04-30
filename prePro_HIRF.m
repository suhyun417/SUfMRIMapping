% clearing variables in the work space
clear all;
close all;

mkdir('figures');

% define file names to read
% expanding ring runs
% fName_R{1}='1_EXP1.nii';
% fName_R{2}='5_EXP2.nii';
% fName_Save='data_expR_mg';

% contracting ring runs
% fName_R{1}='2_CONT1.nii';
% fName_R{2}='6_CONT2.nii';
% fName_Save='data_contR_mg';

% cw-rotating wedge
% fName_R{1}='3_CW1.nii';
% fName_R{2}='7_CW2.nii';
% fName_Save='data_cwWdg_mg';

% ccw-rotating wedge
% fName_R{1}='4_CCW1.nii';
% fName_R{2}='8_CCW2.nii';
% fName_Save='data_ccwWdg_mg';

% HIRF runs
fName_R{1}='9_HRF.nii';
fName_Save='data_HIRF';


d_R=[];
h_R=[];

ps_R=[];
ps_SD=[];

for iRun=1:length(fName_R)
    [d_R{iRun},h_R{iRun}]=cbiReadNifti(fName_R{iRun},{[],[],[],[]},'native',1);
end
% options for high-pass filtering (getting rid of slow drfit fo MR signals)
optDetrend = 1; %1, default/recommended (cutoff frequency=1/128 Hz); 2, 1/(2*cycle duration) Hz for block/periodic design; 

% building up a temporal filter for detreding raw time series
[nR,nC,nS,nT]=size(d_R{1});  % nR ...= # of Rows, Columns, Slices, Time frames
durFrame_sec=1.5; % sampling unit in time for fMRI measurements

numCycle=8; % # of cycles (block design)
numCycle_discard=0; % # of cycles that should be discarded
numCycle_valid=numCycle-numCycle_discard;
numFramePerCycle=nT/numCycle;
numFrame_valid=numCycle_valid*numFramePerCycle;
        
switch (optDetrend)
    case 1,
        FL_frame=nT/(300/durFrame_sec); %nT/(128/durFrame_sec); % default (cutoff frequency=1/128 Hz, conventional)
    case 2,
        FL_frame=nT/((2*numFramePerCycle)/durFrame_sec); % 1/(2*cycle duration) Hz for block/periodic design
end

% Making a Butterworth filter (a.k.a. "maximally flat magnitude filter") for high-pass filtering (which will be used to get rid of
% low-temporal frequency fluctuations) 
% see http://en.wikipedia.org/wiki/Butterworth_filter for more information
% about the BW filter
N=8; % order of filter
xx_radian=(2*pi)/numFrame_valid:(2*pi)/numFrame_valid:(2*pi); 
xx_sec  = [durFrame_sec:durFrame_sec:durFrame_sec*nT]-durFrame_sec/2;
xx_frame		=	[0:1:numFrame_valid*2-1]-numFrame_valid; % frame axis
xx_Hz       = (xx_frame/nT)*(1/(2*durFrame_sec)); 
FL_Hz=(FL_frame/nT)/durFrame_sec;

BWFhp=1./sqrt(1+(2*FL_frame./xx_frame).^(2*N));


% Demonstration of filtering effects
FG1=figure(111); clf; 
SP=subplot(3,3,1), hold on;
title('Low pass filter');
plot(xx_Hz,BWFhp,'k.-','LineWidth',1, 'MarkerSize',6);
xlim([0 max(xx_Hz)]); xlabel('Temporal frequency (Hz)')
ylim([0 1.2]); ylabel('Gain')
set(line([0 0]+FL_Hz, [0 1.2]),'Color','m', 'LineWidth',1)
set(text(FL_Hz+.02, .5, ['Cutoff frequency = ', num2str(FL_Hz), ' Hz; Filter order = ', num2str(N)]), 'Color','m')

SP=subplot(3,3,2),hold on;
title('True time series: 8cycle sine function + synthetic white noise with DC=100');
valDC=100;
snr=10;
ts_true=snr.*sin(numCycle.*xx_radian')+randn(nT,1)+valDC; 
plot(xx_sec, ts_true,'b.-','LineWidth',1, 'MarkerSize',6);
xlim([min(xx_sec) max(xx_sec)]); xlabel('Time (S)')
ylim([-1.5 1.5]*snr+valDC); ylabel('Measurement (a.u.)');

SP=subplot(3,3,3),hold on;
title('Drift time series: synthetic cosine function');
sdr=3;
ts_drift=sdr.*(2*cos(.2*xx_radian')+(2/3)*cos(.6*xx_radian')+(2/5)*cos(1.0*xx_radian'));
ts_drift=ts_drift-mean(ts_drift);
plot(xx_sec, ts_drift,'r-','LineWidth',1);
xlim([min(xx_sec) max(xx_sec)]); xlabel('Time (S)')
ylim([-1.5 1.5]*snr); ylabel('Measurement (a.u.)');

SP=subplot(3,3,4),hold on;
title('Raw time series: synthetic white noise + synthetic cosine function');
ts_raw=ts_true+ts_drift;
plot(xx_sec, ts_raw,'k.-','LineWidth',1, 'MarkerSize',6);
xlim([min(xx_sec) max(xx_sec)]); xlabel('Time (S)')
ylim([-1.5 1.5]*snr+valDC); ylabel('Measurement (a.u.)');

ts=[];
ts(1:numFrame_valid/2,1)=linspace(mean(ts_raw),ts_raw(1), numFrame_valid/2);%ts_org(numImage_valid/2+1:numImage_valid);
ts(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid,1)=ts_raw;
ts(numFrame_valid/2+numFrame_valid+1:numFrame_valid*2,1)=linspace(ts_raw(end), mean(ts_raw), numFrame_valid/2)%ts_org(1:numImage_valid/2);

% lowpass filtering
ts_norm=ts-mean(ts);
Fd		=	fft(double(ts_norm'));
Ms		=	abs(Fd);
Ps		=	angle(Fd);

SP=subplot(3,3,5), cla, hold on;
title('FFT of raw time series');
plot(xx_Hz,fftshift(Ms),'k.-','LineWidth',1, 'MarkerSize',6);
xlim([0 max(xx_Hz)]); xlabel('Temporal frequency (Hz)')
ylabel('Amplitude')
set(line([0 0]+FL_Hz, [0 max(Ms)]),'Color','m', 'LineWidth',1)


BWFbpShift=ifftshift(BWFhp);
NewMs	=	Ms.*BWFbpShift;

SP=subplot(3,3,6), cla, hold on;
title('FFT filtered');
plot(xx_Hz,fftshift(NewMs),'m.-','LineWidth',1, 'MarkerSize',6);
xlim([0 max(xx_Hz)]); xlabel('Temporal frequency (Hz)')
ylabel('Amplitude')
set(line([0 0]+FL_Hz, [0 max(Ms)]),'Color','k', 'LineWidth',1)

compvect=abs(NewMs).*( cos(Ps)+j.*sin(Ps) );
vect=ifft(compvect);
BpData=real(vect);

ts_new_padded=mean(ts)+BpData;

% convert into percent change time series
ps_new_padded=100*(ts_new_padded-mean(ts_new_padded))/mean(ts_new_padded);
ps_new=ps_new_padded(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid);

ps_true=100*( (ts_true-mean(ts_true))./mean(ts_true));
ps_raw=100*( (ts_raw-mean(ts_raw))./mean(ts_raw));

            
SP=subplot(3,3,7), cla; hold on;
title('Filtered time series');
plot(xx_sec, ps_new,'m.-','LineWidth',1, 'MarkerSize',6);
plot(xx_sec, ps_raw-ps_new','r.-','LineWidth',1, 'MarkerSize',6);
xlim([min(xx_sec) max(xx_sec)]); xlabel('Time (S)')
ylim([-1.5 1.5]*snr); ylabel('Measurement (a.u.)');




% comparison between True versus Raw
SP=subplot(3,3,8), cla; hold on;
title('True versus Raw');

plot(ps_true, ps_raw, 'k.', 'MarkerSize', 6);
axis equal; axis square;
set(line([-1 1]*max([max(abs(ps_true)),max(abs(ps_raw))]), [-1 1]*max([max(abs(ps_true)),max(abs(ps_raw))])), 'LineStyle','--','Color', 'k');
xlim([-1 1]*max([max(abs(ps_true)),max(abs(ps_raw))])); xlabel('True TS')
ylim([-1 1]*max([max(abs(ps_true)),max(abs(ps_raw))])); ylabel('Raw TS')

% comparison between True versus Filtered
SP=subplot(3,3,9), cla; hold on;
title('True versus Raw');
plot(ps_true, ps_new, 'm.', 'MarkerSize', 6);
axis equal; axis square;
set(line([-1 1]*max([max(abs(ps_true)),max(abs(ps_raw))]), [-1 1]*max([max(abs(ps_true)),max(abs(ps_raw))])), 'LineStyle','--','Color', 'k');
xlim([-1 1]*max([max(abs(ps_true)),max(abs(ps_raw))])); xlabel('True TS')
ylim([-1 1]*max([max(abs(ps_true)),max(abs(ps_raw))])); ylabel('Filtered TS')


% Tip: you can modfy the cut frequency & smoothness of the filter by
% varying FL and N

% definition of histogram parameters
numHist=100;
cut_yHist=200;

FG2=figure(2); clf;
set(FG2,'Name', 'Mean Intensity Images of Functionals');

FG3=figure(3); clf;
set(FG3,'Name', 'Distribution of Mean Intensity Values');

FG4=figure(4); clf;
set(FG4,'Name', 'Pooled Distribution of Mean Intensity Values');

for iRun=1:length(fName_R)
    % show the mean EPI averaged across time

    for iS=1:nS
        meanEPI_acTime{iRun}(:,:,iS)=mean(d_R{iRun}(:,:,iS,:),4);
    end
    maxVal=max(max(max(meanEPI_acTime{iRun}))); % for divisive normalization for image intensities
    minVal=1;



%     xHist=[(maxVal/numHist):(maxVal/numHist):maxVal]-(maxVal/numHist)/2;
    xHist_log10=linspace(log10(minVal), log10(maxVal), numHist);
    yHist_stack{iRun}=[];

    countSlice=0;
    rotVal=-1;%-1; % rotation parameter for image illustration
    tic
    for i=1:4
        for ii=1:6
            countSlice=countSlice+1;
            fprintf(1,'Slice:%d elapsed time: %3.1f\n',  countSlice, toc);

            if( iRun==1)
                figure(FG2);
                SP=subplot(4,6,countSlice); hold on;
                image(255*(rot90(meanEPI_acTime{iRun}(:,:,countSlice)./maxVal,rotVal)));
                colormap(gray(256));
                axis equal; 
                title(['Slice #',num2str(countSlice)]);

                [tempR,tempC]=size(255*(rot90(meanEPI_acTime{iRun}(:,:,countSlice)./maxVal,rotVal)));
                xlim([1 tempC]); xlabel('Voxel')
                ylim([1 tempR]); ylabel('Voxel')
            end

            figure(FG3);
            SP=subplot(4,6,countSlice); hold on;
            title(['Slice #',num2str(countSlice)]); hold on;
%             yHist=hist(reshape(meanEPI_acTime{iRun}(:,:,countSlice),nC*nR,1),xHist);
            yHist_log10=hist(log10(reshape(meanEPI_acTime{iRun}(:,:,countSlice),nC*nR,1)),xHist_log10);
            yHist_stack{iRun}(countSlice,:)=yHist_log10;

            PL=plot(xHist_log10, yHist_log10,'k-');
            if(iRun==2) set(PL, 'LineStyle','--'); end; 
            xlim([1 max(xHist_log10)]); xlabel('Image intensity (log10)');
            ylim([0 cut_yHist]); ylabel('Number of voxels');
            drawnow;
        end
    end

    figure(FG4);
    SP=subplot(2,2,1); hold on;
    set(FG4,'Name', 'Summed Distribution of Mean Intensity Values');
    PL=plot(xHist_log10, sum(yHist_stack{iRun}),'k-');
    if(iRun==2) set(PL, 'LineStyle','--'); end; 
    xlim([1 max(xHist_log10)]); xlabel('Image intensity (log10)');
    ylim([0 cut_yHist*24]); ylabel('Number of voxels');
    drawnow;

    SP=subplot(2,2,2); hold on;
    cmap_slice=jet(nS);
    for iS=1:nS
        yProp=yHist_stack{iRun}(iS,:)/sum(yHist_stack{iRun}(iS,:));
        if(iS>1 & iS <19)
            PL=plot(xHist_log10, yProp,'-','Color',cmap_slice(iS,:),'LineWidth',1); hold on;
             if(iRun==2) set(PL, 'LineStyle','--'); end;
        end
        vectIntAvg(iS,1)=sum(yProp.*xHist_log10);
    end
    xlim([1 max(xHist_log10)]); xlabel('Image intensity (log10)');
    ylim([0 .03]); ylabel('Proportion of voxels');
    drawnow;


    % filtering time series & convesion to Percent Change Signals
    % (normalization for inhomogeniety in field strength, tissue, local magnetic susceptibility, distance from the receiver coil etc )
    tic
    for iS=1:nS
        fprintf(1,'Slice:%d elapsed time: %3.1f\n',iS, toc);
        for iR=1:nR

            for iC=1:nC

            ts_org=squeeze(d_R{iRun}(iR,iC,iS,numCycle_discard*numFramePerCycle+1:nT));

            % bufferring the head and tail part of the original time series (to
            % get rid of artifactual temporal frequency components due to
            % temporal edges)
            ts(1:numFrame_valid/2,1)=mean(ts_org);%ts_org(numImage_valid/2+1:numImage_valid);
            ts(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid,1)=ts_org;
            ts(numFrame_valid/2+numFrame_valid+1:numFrame_valid*2,1)=mean(ts_org);%ts_org(1:numImage_valid/2);

            % lowpass filtering
            Fd		=	fft(double(ts'));
            Ms		=	abs(Fd);
            Ps		=	angle(Fd);

            BWFbpShift=ifftshift(BWFhp);
            NewMs	=	Ms.*BWFbpShift;
            
            figure;
            title('FFT of raw time series');
            plot(xx_Hz,fftshift(Ms),'k.-','LineWidth',1, 'MarkerSize',6);
            xlim([0 max(xx_Hz)]); xlabel('Temporal frequency (Hz)')
            ylabel('Amplitude')
            set(line([0 0]+FL_Hz, [0 max(Ms)]),'Color','m', 'LineWidth',1)

            compvect=abs(NewMs).*( cos(Ps)+j.*sin(Ps) );
            vect=ifft(compvect);
            BpData=real(vect);

            ts_new=mean(ts)+BpData;

            % convert into percent change time series
            ps_new=100*(ts_new-mean(ts_new))/mean(ts_new);

            ps_R{iRun}(iR,iC,iS,:)=ps_new(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid);
            ps_SD{iRun}(iR,iC,iS)=std(ps_new(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid));
            
            estDrift=double(100*(ts_org'-mean(ts_org))/mean(ts_org))-ps_new(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid);
            [P,S]=polyfit(1:nT,estDrift,1);
            ps_DriftSlope{iRun}(iR,iC,iS)=P(1);
    %         plot(ts_new(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid),'r-'); drawnow; pause;
    %         plot([100*(ts-mean(ts_new))/mean(ts_new),ps_new']); drawnow; pause;
            end
        end
    end
    
    SP=subplot(2,2,3); hold on;
    
    xHist_ds=-.2:.005:.2;
    for iS=1:nS
        vectDS= reshape(ps_DriftSlope{iRun}(:,:,iS), nR*nC,1);
        yHist_DS=hist(vectDS, xHist_ds); 
        PL=plot(xHist_ds, yHist_DS,'-', 'MarkerSize', 6, 'Color',cmap_slice(iS,:),'LineWidth',1); hold on;
        if(iRun==2) set(PL, 'LineStyle','--'); end; 
    end
    xlim([min(xHist_ds) max(xHist_ds)]); xlabel('Drift slope (a.u.)');
    ylabel('Number of voxels');
    drawnow;
    
    SP=subplot(2,2,4); hold on;
    
    for iS=1:nS
        vectDS= reshape(ps_DriftSlope{iRun}(:,:,iS), nR*nC,1);
        vectSD= reshape(ps_SD{iRun}(:,:,iS), nR*nC,1);
        PL=plot(vectDS, log10(vectSD),'o', 'MarkerSize', 4, 'Color',[0 0 0],'MarkerFaceColor', [0 0 0]); hold on;
        if(iRun==2) set(PL,'MarkerFaceColor',[0 0 0]+.5); end
    end
    xlim([-1 1]*1.5); xlabel('Drift slope (%signal/frame)');
    ylim([0 2.5]); ylabel('SD in % signal (log10)');
    drawnow;
    
end


% Identify blood-clamping voxels by applying the "SD>10%" rule
% show distributions of SD values
FG7=figure(7); clf;
set(FG7,'Name', 'Distrbution of SD (%)');
maxSD_log=2;
minSD_log=-.5;
xHist=linspace(-.5,maxSD_log,200);

global obsData flagDraw
flagDraw=0;
obsData=[];
guessIn = [0.4 1.12 .3 .05 .27]; %mu1 mu2 sig1 sig2 ratio1
pVal_critBlood=10^-8;

for iRun=1:length(fName_R)
    SP=subplot(2,2,iRun); cla; hold on;
    if(iRun==2) guessIn=guessOut; end;
     %max([max(max(max(ps_SD{1}))), max(max(max(ps_SD{2})))]);
    
    yHist=hist(reshape(log10(ps_SD{iRun}), nR*nC*nS,1), xHist);
    plot(xHist,yHist./sum(yHist),'k-');
    xlim([minSD_log maxSD_log]); xlabel('log10 of SD (% singal)');
    ylim([0 .1]); ylabel('Proportion');
    
    % two gaussian fitting
    vectLog10SD_sort=sort(reshape(log10(ps_SD{iRun}), nR*nC*nS,1));
    vectTemp=vectLog10SD_sort(vectLog10SD_sort>minSD_log & vectLog10SD_sort<maxSD_log);
    obsData=[vectTemp, [1:1:length(vectTemp)]'./length(vectTemp)];
    
    tic
    guessOut=fminsearch('sse2cdf', guessIn)
    toc
    matParam(iRun,:)=guessOut;
    
    mu1=guessOut(1);
    mu2=guessOut(2);
    sig1=guessOut(3);
    sig2=guessOut(4);
    ratio1=guessOut(5);
% 
%     y1=ratio1*(normpdf(obsData(:,1),mu1,sig1)./normpdf(0,mu1,sig1));
%     y2=(1-ratio1)*(normpdf(obsData(:,1),mu2,sig2)./normpdf(0,mu2,sig2));
    
     y1=ratio1*normpdf(xHist,mu1,sig1);
    y2=(1-ratio1)*normpdf(xHist,mu2,sig2);
    y11=y1./sum(y1+y2);
    y22=y2./sum(y1+y2);
    figure(FG7)
    plot(xHist,y11,'b-','LineWidth',2);
    plot(xHist,y22,'r-','LineWidth',2);
    plot(xHist,yHist./sum(yHist),'k-'); drawnow;
end
% define criterion
MuSig2_avg=mean(matParam(:,[2 4]),1);
critBlood=norminv(pVal_critBlood, MuSig2_avg(1), MuSig2_avg(2) );
for iRun=1:length(fName_R)
    SP=subplot(2,2,iRun);
    set(line([0 0]+critBlood, [0 .1]),'LineWidth',2, 'LineStyle', '--', 'Color', 'm');
end


SP=subplot(2,2,3); cla; hold on;

temp1=reshape(log10(ps_SD{1}), nR*nC*nS,1);
tLoc1=find(temp1>=critBlood);

plot(temp1,temp1,'bo');
plot(temp1(tLoc1),temp1(tLoc1),'ro');
xlim([minSD_log maxSD_log]); xlabel('SD (% singal): run1');
ylim([minSD_log maxSD_log]); ylabel('SD (% singal): run1');

% plotting SD values as a function of mean intensity values
SP=subplot(2,2,4); cla; hold on;
vectSD_pooled= reshape(log10(sqrt( (ps_SD{1}.^2+ps_SD{1}.^2)/2)), nR*nC*nS,1);
vectImgInt_pooled=reshape((meanEPI_acTime{1}+meanEPI_acTime{1})/2, nR*nC*nS,1);
plot(log10(vectImgInt_pooled), vectSD_pooled,'r.','MarkerSize',2);
plot(log10(vectImgInt_pooled(tLoc1)), vectSD_pooled(tLoc1),'b.','MarkerSize',2);

set(line([1 3.5], [0 0]+critBlood),'LineWidth',2, 'LineStyle', '--', 'Color', 'm');
xlim([1 3.5]); xlabel('Pooled image intensity (log10)');
ylim([minSD_log 2]); ylabel('log10 of Pooled SD (% singal)');

% Making a valid brain region mask
mask3D_BloodProof=log10(ps_SD{1})<critBlood & log10(ps_SD{1})<critBlood;
countSlice=0;
rotVal=1; % rotation parameter for image illustration
FG8=figure(8); clf;
meanEPI_pooled=(meanEPI_acTime{1}+meanEPI_acTime{1})/2;

for i=1:4
    for ii=1:6
        countSlice=countSlice+1;
        fprintf(1,'Slice:%d elapsed time: %3.1f\n',  countSlice, toc);

%         image(255*(meanEPI_acTime(:,:,countSlice)./maxVal));
        figure(FG8);
        SP=subplot(4,6,countSlice); 
        temp=meanEPI_pooled(:,:,countSlice).*mask3D_BloodProof(:,:,countSlice);
        image(255*(rot90(temp./maxVal,rotVal)));
        colormap(gray(256));
        axis equal; 
        title(['Slice #',num2str(countSlice)]);

        [tempR,tempC]=size(255*(rot90(temp)));
        xlim([1 tempC]); xlabel('Voxel')
        ylim([1 tempR]); ylabel('Voxel')
    end
end



raw_mg= d_R{1};
ps_mg= ps_R{1};
sd_mg= ps_SD{1};
            

save(fName_Save, 'raw_mg', 'ps_mg', 'sd_mg', 'mask3D_BloodProof', 'h_R')



%%
saveas(FG1, ['figures/' mfilename '_figure1.fig']);
saveas(FG2, ['figures/' mfilename '_figure2.fig']);
saveas(FG3, ['figures/' mfilename '_figure3.fig']);
saveas(FG4, ['figures/' mfilename '_figure4.fig']);
saveas(FG7, ['figures/' mfilename '_figure7.fig']);
saveas(FG8, ['figures/' mfilename '_figure8.fig']);

