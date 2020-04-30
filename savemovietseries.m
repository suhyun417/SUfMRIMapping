% savemovietseries
% save movie tseries of fMRI data

function savemovietseries()
  
  addpath ('/projects/parksh') %addpath ('/einstein0/USRlab/projects/parksh')

%   global S dataBOLD 
%   global STDPATH DSP DATA GH   
%     
%   % Load a data file containing neural regressor and fMRI time series,
%   % along with correlation maps for each movie
%   
%   % First start,open the XL file, and reset the windows
%   blockpaths_soohyun;
%   STDPATH.xlsfile = 'block_timecourse_files_ber.xls';
%   loadXLSFile;
% %  resetBlockanaWindows;
%   
%   % File params
%   remakeMaps = 0;   
%   nameSubjNeural = 'Tor'; %'Rho'; % Monkey name for neural regressor
  nameSubjBOLD = 'Art'; %'Ava'; % Monkey name for fMRI data
%   STDPATH.dataNeural = fullfile(STDPATH.datahome, nameSubjNeural);
  dirDataBOLD = fullfile(STDPATH.datahome, nameSubjBOLD);
%   
%   saveFileName_neural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
  saveFileName_BOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

%   saveFileName = ['dataNeuralBOLD_', nameSubjNeural, nameSubjBOLD, '_indMov.mat']; %'dataNeuralBOLD_TorArt_3Movies'; %'dataNeuralBOLD_TorArt_123'; %'dataNeuralBOLD_TorArt_indMov'; %'dataNeuralBOLD_TorArt'; %'dataNeuralBOLD_TorArt_131415'; %'dataNeuralBOLD_TorArt'; %'dataNeuralBold';
%   saveFileName_matSDF = [nameSubjNeural, '_tempMatSDF_indMov.mat'];
  
  skip=7; % number of TRs to skip 
  redo=0;
  textfileroot = 'FL_e66_9FILES';  %% MUST be matched to the fMRI data
 
  MCD = loadMeanConcatData(textfileroot,skip,redo);
  
  % save file to directory
  save(fullfile(dirDataBOLD, saveFileName_BOLD), 'MCD') 

  

function MCD = loadMeanConcatData(textfileroot,skip,redo)
  
  global STDPATH
  
  sessdatafile = sprintf('%s/%s.mat',STDPATH.sessdatdir,textfileroot);
  sesstextfile = sprintf('%s.txt',textfileroot);
  SI = SU_createSessInfo(sesstextfile,[],skip); % get infor from text file
  filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);
  % 
  % Init by loading one file
  %
  s_sub = filelist{1}.subj;
  s_ses = filelist{1}.sess;
  s_sc  = filelist{1}.scan;
   [fmri_tc_3d,dgz] = SU_loadFile(s_sub,s_ses,s_sc);
  
  [TR,NR] = getTRandNR();
  
  if ~isempty(dir(sessdatafile)) & ~redo
	load(sessdatafile);  % load already created file
  else
	% Create and save the mean, concatenated files
	%
	MCD = buildMeanConcatData(SI);     % create the big S structure
 	save(sessdatafile,'MCD');            % save the file
  end
  
  
function MCD = buildMeanConcatData(SI)
  
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
    
	for f=1:totfiles
	  s_sub = filelist{f}.subj;
	  s_ses = filelist{f}.sess;
	  s_sc  = filelist{f}.scan;
	  
	  [fmri_tc,dgz,notes] = SU_loadFile(s_sub,s_ses,s_sc);

	  
	  if f==1
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
	  allvoltc(:,:,:,:,f) = clp_fmri_tc;
	  alldgz(f) = dgz;
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
%	dt_mvoltc = SU_detrendMovDat(mvoltc);
     dt_mvoltc = mvoltc;
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
  
