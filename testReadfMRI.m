% testReadfMRI
% savetseries

addpath ('/einstein0/USRlab/projects/parksh/analysis')

global S dataBOLD
global STDPATH DSP DATA GH

% First start,open the XL file, and reset the windows
blockpaths_soohyun;
STDPATH.xlsfile = 'block_timecourse_files_ber.xls';
loadXLSFile;

textfileroot = 'FL_e66_9FILES'; %'FL_e66_allfiles2'; %% MUST be matched to the fMRI data

%   MCD = loadMeanConcatData(textfileroot,skip,redo);

%   sessdatafile = sprintf('%s/%s.mat',STDPATH.sessdatdir,textfileroot);
sesstextfile = sprintf('%s.txt',textfileroot);
SI = SU_createSessInfo(sesstextfile,[],skip); % get infor from text file
filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);

% Init by loading one file
s_sub = filelist{1}.subj;
s_ses = filelist{1}.sess;
s_sc  = filelist{1}.scan;
[fmri_tc_3d,dgz] = SU_loadFile(s_sub,s_ses,s_sc);

d_R{iRun} = fmri_tc_3d;

for iS=1:nS
    
    for iR=1:nR
        
        for iC=1:nC
            fprintf(1,'Slice:%d/%d Column:%d/%d Row:%d/%d \n',iS, nS, iR, nR, iC, nC);
            
            ts_org=squeeze(d_R{iRun}(iR,iC,iS,nT-numFrame_valid+1:nT)); %squeeze(d_R{iRun}(iR,iC,iS,numCycle_discard*numFramePerCycle+1:nT));
            
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
            
            figure(105); clf;
            title('FFT of raw time series');
            plot(xx_Hz,fftshift(Ms),'k.-','LineWidth',1, 'MarkerSize',6);
            xlim([0 max(xx_Hz)]); xlabel('Temporal frequency (Hz)')
            ylabel('Amplitude')
            set(line([0 0]+FL_Hz, [0 max(Ms)]),'Color','m', 'LineWidth',1)
            input('')
            
        end
    end
end




