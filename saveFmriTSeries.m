% saveFmriTSeries: nameSubj, optDetrend, optPercentSig
% save two types of files
%  --1) file per each movie, containing each run
%  --2) one file containing averaged across runs for each movie


addpath ('/einstein0/USRlab/projects/parksh/analysis')

%% File Details
nameSubj = 'Art'; %'Ava'; %'Art'; % one of the inputs of this function
monkID = SU_getMonkIDByName(nameSubj);

% set save file names
saveFileNameHeader = sprintf('%s_DataFMRI', nameSubj);
if optDetrend
    saveFileNameHeader = [saveFileNameHeader '_det'];
end
if optPercentSig
    saveFileNameHeader = [saveFileNameHeader '_ps'];
end

global S dataBOLD
global STDPATH DSP DATA GH

% First start,open the XL file, and reset the windows
blockpaths_soohyun;
STDPATH.xlsfile = 'block_timecourse_files_ber.xls';
loadXLSFile;


textfileroot = sprintf('FL_%s_allfiles2', monkID);
% textfileroot = 'FL_e66_allfiles2'; %'FL_e66_9FILES'; %'FL_e66_allfiles2'; %% MUST be matched to the fMRI data

%   MCD = loadMeanConcatData(textfileroot,skip,redo);

skip=0;
sesstextfile = sprintf('%s.txt',textfileroot);
SI = SU_createSessInfo(sesstextfile,[],skip); % get infor from text file
filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);


%% For each movie, load each run and detrend it (just for the linear trend in , then change it into percent signal
% Parameters for filtering (detrending)
durFrame_sec=2.4; %1.5; % sampling unit in time for fMRI measurements

% numCycle=8; % # of cycles (block design)
% numCycle_discard=0; % # of cycles that should be discarded
% numCycle_valid=numCycle-numCycle_discard;
% numFramePerCycle=nT/numCycle;
numFrame_valid=nT-7; %numCycle_valid*numFramePerCycle;

FL_frame=1.2; %nT/(300/durFrame_sec); %nT/(128/durFrame_sec); % default (cutoff frequency=1/128 Hz, conventional)

% Making a Butterworth filter (a.k.a. "maximally flat magnitude filter") for high-pass filtering (which will be used to get rid of
% low-temporal frequency fluctuations) 
% see http://en.wikipedia.org/wiki/Butterworth_filter for more information
% about the BW filter
N=8; % order of filter
xx_radian=(2*pi)/numFrame_valid:(2*pi)/numFrame_valid:(2*pi); 
xx_sec  = [durFrame_sec:durFrame_sec:durFrame_sec*numFrame_valid]-durFrame_sec/2;
xx_frame		=	[0:1:numFrame_valid*2-1]-numFrame_valid; % frame axis
xx_Hz       = (xx_frame/nT)*(1/(2*durFrame_sec)); 

FL_Hz=(FL_frame/nT)/durFrame_sec;
BWFhp=1./sqrt(1+(2*FL_frame./xx_frame).^(2*N));

%
unimov = sort(unique(SI.movID));
nunimov = length(unimov);
catimgdat = [];
rgr = zeros(1,nunimov*MAX_TR);
procmovs = [];

for u = 1:nunimov
    movindxs = find(SI.movID == unimov(u));
    nmovindxs = length(movindxs);
    
    % collect all the data files for this movie
    %
    filelist = SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs),...
        SI.scanID(movindxs));
    
    totfiles = length(filelist) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for iFile=1:totfiles
        
        s_sub = filelist{iFile}.subj;
        s_ses = filelist{iFile}.sess;
        s_sc  = filelist{iFile}.scan;
        
        cursession = sprintf('%s.%s',s_sub,s_ses);
        sessindx   = strmatch(cursession,DSP.filedetails.sess);
        scanindx   = find(DSP.filedetails.scan == s_sc);
        xlsindx    = intersect(sessindx,scanindx);
        
        DSP.filedetails.date{xlsindx}=ri;h
        
        [fmri_tc,dgz,notes] = SU_loadFile(s_sub,s_ses,s_sc);
        
        %         % Init by loading one file
        %         s_sub = filelist{1}.subj;
        %         s_ses = filelist{1}.sess;
        %         s_sc  = filelist{1}.scan;
        %         [fmri_tc_3d,dgz] = SU_loadFile(s_sub,s_ses,s_sc);
        
%         d_R{iFile} = fmri_tc_3d;
        [nR,nC,nS,nT]=size(fmri_tc_3d); %size(d_R{iFile});
        
        for iS=1:nS
            
            for iR=1:nR
                
                for iC=1:nC
                    fprintf(1,'Slice:%d/%d Column:%d/%d Row:%d/%d \n',iS, nS, iR, nR, iC, nC);
                    
                    ts_org=squeeze(fmri_tc_3d(iR,iC,iS,nT-numFrame_valid+1:nT)); %squeeze(d_R{iRun}(iR,iC,iS,numCycle_discard*numFramePerCycle+1:nT));
                    
                    % bufferring the head and tail part of the original time series (to
                    % get rid of artifactual temporal frequency components due to
                    % temporal edges)
                    ts(1:numFrame_valid/2,1)=mean(ts_org);%ts_org(numImage_valid/2+1:numImage_valid);
                    ts(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid,1)=ts_org;
                    ts(numFrame_valid/2+numFrame_valid+1:numFrame_valid*2,1)=mean(ts_org);%ts_org(1:numImage_valid/2);
                    
                    % high-pass filtering
                    Fd		=	fft(double(ts'));
                    Ms		=	abs(Fd);
                    Ps		=	angle(Fd);
                    
                    BWFbpShift=ifftshift(BWFhp);
                    NewMs	=	Ms.*BWFbpShift;                    
                    
                    compvect=abs(NewMs).*( cos(Ps)+j.*sin(Ps) );
                    vect=ifft(compvect);
                    BpData=real(vect);
                    
                    ts_new_padded=mean(ts)+BpData;
                    ts_new_bwfilter = ts_new_padded(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid);
                    
                    % convert into percent change time series
                    ps_new_padded=100*(ts_new_padded-mean(ts_new_padded))/mean(ts_new_padded);
                    ps_new_bwfilter=ps_new_padded(numFrame_valid/2+1:numFrame_valid/2+numFrame_valid);
                    
                    %             ps_true=100*( (ts_true-mean(ts_true))./mean(ts_true));
                    ps_org=100*( (ts_org-mean(ts_org))./mean(ts_org));
                    
                    % Linear detrend
                    ts_new_linDetrend = detrend(ts_org);
                    ps_new_linDetrend = 100*(ts_new_linDetrend-mean(ts_new_linDetrend))/mean(ts_new_linDetrend);
                    
                    d_linDetrend_ps{iFile}(iR, iC, iS, :) = ps_new_linDetrend;
                    d_bwfilter_ps{iFile}(iR, iC, iS, :) = ps_new_bwfilter;
                    d_org{iFile}(iR, iC, iS, :) = ps_org;
                    
                    
%                     figure(105); clf;
%                     subplot(3,1,1), cla, hold on;
%                     title('FFT of raw time series');
%                     plot(xx_Hz,fftshift(Ms),'k.-','LineWidth',1, 'MarkerSize',6);
%                     xlim([0 max(xx_Hz)]); xlabel('Temporal frequency (Hz)')
%                     ylabel('Amplitude')
%                     set(line([0 0]+FL_Hz, [0 max(Ms)]),'Color','m', 'LineWidth',1)
%                     
%                     subplot(3,1,2), cla, hold on;
%                     title('Raw and filtered time series: using BW filter');
%                     plot(xx_sec, ts_org,'b.-','LineWidth',1, 'MarkerSize',6);
%                     plot(xx_sec, ts_new_bwfilter,'m.-','LineWidth',1, 'MarkerSize',6);
%                     plot(xx_sec, ts_org-ts_new_bwfilter'+mean(ts),'r.-','LineWidth',1, 'MarkerSize',6);
%                     xlim([min(xx_sec) max(xx_sec)]); xlabel('Time (s)')
%                     ylabel('Measurement (a.u.)')
%                     
%                     subplot(3,1,3), cla; hold on;
%                     title('Raw and filtered time series: linear detrending');
%                     plot(xx_sec, ts_org,'b.-','LineWidth',1, 'MarkerSize',6);
%                     plot(xx_sec, ts_new_linDetrend+mean(ts),'m.-','LineWidth',1, 'MarkerSize',6);
%                     plot(xx_sec, ts_org-ts_new_linDetrend,'r.-','LineWidth',1, 'MarkerSize',6);
%                     xlim([min(xx_sec) max(xx_sec)]); xlabel('Time (s)')
%                     ylabel('Measurement (a.u.)');
%                     %             input('')
%                     
%                     
%                     %             subplot(3,1,3), cla; hold on;
%                     %             title('Filtered time series in percent signal');
%                     %             plot(xx_sec, ps_new,'m.-','LineWidth',1, 'MarkerSize',6);
%                     %             plot(xx_sec, ps_raw-ps_new','r.-','LineWidth',1, 'MarkerSize',6);
%                     %             xlim([min(xx_sec) max(xx_sec)]); xlabel('Time (s)')
%                     %             ylabel('Measurement (a.u.)');
                    
                end
            end
        end
        % save data to a separate file for each movie containing all individual
        % runs
        
    end
    % Average data across runs for each movie   
        
    % Put it into a struct
    
   
end
% Save data averaged across runs into another file




%%
MAX_TR = SI.max_tr;
  skip=SI.skiptr;
  
  unimov = sort(unique(SI.movID));
  nunimov = length(unimov);
  catimgdat = [];
  rgr = zeros(1,nunimov*MAX_TR);
  procmovs = [];

  for u = 1:nunimov
	movindxs = find(SI.movID == unimov(u));
	nmovindxs = length(movindxs);
	
	% collect all the data files for this movie
	%
	filelist = SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs),...
							   SI.scanID(movindxs));

	totfiles = length(filelist);
    
	for iFile=1:totfiles
	  s_sub = filelist{iFile}.subj;
	  s_ses = filelist{iFile}.sess;
	  s_sc  = filelist{iFile}.scan;
	  
	  [fmri_tc,dgz,notes] = SU_loadFile(s_sub,s_ses,s_sc);
	  
	  if iFile==1
		[a,b,c,t] = size(fmri_tc);
		allvoltc = zeros(a,b,c,MAX_TR,totfiles);
	  end
	  
	  scanlen = size(fmri_tc,4);
	  if scanlen == 250
		mov_start_tr = 63;% offset when movie starts in scan
	  elseif scanlen == 125
		mov_start_tr = 1; % short movie
	  else
		fprintf('WARNING: bad scan length %d\n',scanlen);
      end
      
	  val_tr         = [mov_start_tr:(mov_start_tr+MAX_TR-1)]; 
	  [clp_fmri_tc,clp_mdgz] = SU_clipMovDat(fmri_tc,dgz,val_tr);
	  allvoltc(:,:,:,:,iFile) = clp_fmri_tc;
	  alldgz(iFile) = dgz;
	end
	[mvoltc,mdgz]=SU_computeMeanTC(allvoltc,alldgz);

    if skip ~= 0
      vox_mns=nanmean(mvoltc,4);
      for i=1:skip
		% sets the skipped TRs to be the mean of the voxel 10/30/12 BER
        % 
		mvoltc(:,:,:,i)=vox_mns(:,:,:); 
      end
    end
    
	fprintf('Starting at vol %d. Detrending mov %d\n',mov_start_tr,unimov(u)); 
	dt_mvoltc = SU_detrendMovDat(mvoltc);

    if skip ~=0
	  % sets the skipped TRs to be NAN 10/30/12 BER 	  
	  %
	  dt_mvoltc(:,:,:,1:skip)=NaN; 
    end
	allmvoltc(:,:,:,:,u) = dt_mvoltc;
	allmdgz(u) = clp_mdgz;
	catimgdat = cat(2,catimgdat,SU_getImageDat(unimov(u)));
  end

  [catmvoltc,catmdgz]  =...
	  SU_concatVolumesAndDGZ(allmvoltc,allmdgz);

  % Store Scan Info in S % Added by BER 02/12/13 
  %
  MCD.max_tr=MAX_TR;
  MCD.Fs=SI.Fs;
  MCD.TR=SI.TR;
  MCD.skiptr=SI.skiptr;
  
  MCD.catmvoltc = catmvoltc;
  MCD.catmdgz   = catmdgz;
  MCD.catimgdat = catimgdat;
  MCD.unimov    = unimov;



