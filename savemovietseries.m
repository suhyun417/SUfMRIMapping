savemovietseries
function S_neuralRegressor()
  
  addpath ('/projects/parksh') %addpath ('/einstein0/USRlab/projects/parksh')

  global S dataBOLD 
  global STDPATH DSP DATA GH   
    
  % Load a data file containing neural regressor and fMRI time series,
  % along with correlation maps for each movie
  
  % First start,open the XL file, and reset the windows
  blockpaths_soohyun;
  STDPATH.xlsfile = 'block_timecourse_files_ber.xls';
  loadXLSFile;
%  resetBlockanaWindows;
  
  % File params
  remakeMaps = 0;   
  nameSubjNeural = 'Tor'; %'Rho'; % Monkey name for neural regressor
  nameSubjBOLD = 'Art'; %'Ava'; % Monkey name for fMRI data
  STDPATH.dataNeural = fullfile(STDPATH.datahome, nameSubjNeural);
  STDPATH.dataBOLD = fullfile(STDPATH.datahome, nameSubjBOLD);
  
  saveFileName_neural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
  saveFileName_BOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

%   saveFileName = ['dataNeuralBOLD_', nameSubjNeural, nameSubjBOLD, '_indMov.mat']; %'dataNeuralBOLD_TorArt_3Movies'; %'dataNeuralBOLD_TorArt_123'; %'dataNeuralBOLD_TorArt_indMov'; %'dataNeuralBOLD_TorArt'; %'dataNeuralBOLD_TorArt_131415'; %'dataNeuralBOLD_TorArt'; %'dataNeuralBold';
%   saveFileName_matSDF = [nameSubjNeural, '_tempMatSDF_indMov.mat'];
  
  skip=7; % number of TRs to skip
  textfileroot = 'FL_e66_9FILES';  %% MUST be matched to the fMRI data
 
%   MCD = loadMeanConcatData(textfileroot,skip,redo);
  
%   sessdatafile = sprintf('%s/%s.mat',STDPATH.sessdatdir,textfileroot);
  sesstextfile = sprintf('%s.txt',textfileroot);
  SI = SU_createSessInfo(sesstextfile,[],skip); % get infor from text file
  filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);
  
%   fullpath = sprintf('%s/%s',STDPATH.filelists, sesstextfile);
%   if ~isempty(dir(fullpath))
% 	fid = fopen(fullpath,'r');
% 	dat = textscan(fid,'%s %s %d','CommentStyle','%');
% 	fclose(fid);
%   else
% 	fprintf('No such file %s\n',fullpath);
% 	return;
%   end
%   
%   subjs = dat{1};
%   sesss = dat{2};
%   scans = dat{3};
  
  % Init by loading one file
  s_sub = filelist{1}.subj;
  s_ses = filelist{1}.sess;
  s_sc  = filelist{1}.scan;
  [fmri_tc_3d,dgz] = SU_loadFile(s_sub,s_ses,s_sc);
  
  if ~remakeMaps
     if exist('S') & exist('dataBOLD') & ~isempty(S) & ~isempty(dataBOLD)
        fprintf('Using S and dataBOLD already in memory\n');
     else
%        fprintf('Loading File %s\n',saveFileName);
%        load(fullfile(STDPATH.datahome,saveFileName), 'S');

       fprintf('Loading File %s and %s\n',saveFileName_neural,saveFileName_BOLD);
       load(fullfile(STDPATH.dataNeural,saveFileName_neural));
       load(fullfile(STDPATH.dataBOLD,saveFileName_BOLD), 'dataBOLD');
     end
  else
    % build the S structure
    [S, dataBOLD] = getCorrMap_indMov(nameSubjNeural, nameSubjBOLD, saveFileName_matSDF, saveFileName); %[S,dataBOLD] = getCorrMap();
  end
  
  % Get movie IDs common in two dataset
  [setMovIDs, locMovBOLD] = intersect(dataBOLD.unimov, cat(2,S(1,:).movID)); 
  dataBOLD.mvoltc = dataBOLD.mvoltc(setMovIDs);
  dataBOLD.unimov = setMovIDs;
  
%   listBestCellInd = [1 5 6 9 11 14 16 21 26 29 31 32 33 35 36 39]; %[1 5 6 9 10 14 15 16 21 23 25 26 31 
%   iChan = 1;% listBestCellInd(2);
%   iMov = 1;
  
%   seeCorrMap(S, dataBOLD, iChan, iMov);
  
%   DSP.proc.scalarmap_3d = S(iChan, iMov).mapR; %S(iChan, iMov).mapR; %.*-1; %S(iChan).mapR.*-1; %S(6).mapR.*-1;
%   DSP.proc.fmri_tc_3d= voltc;

  fillGlobalFieldsNeural(S, dataBOLD) %, iChan); %fillGlobalFieldsNeural(S,MCD,dataBOLD,iChan);

      
  showVolumeTimeCourse;
  setLightboxOrientation('hor');
  updateLightbox;

return

  % setup advanced functional
  set(GH.main.but.advanced,'Value',1);
  doAdvancedFunctional;
  
  set(GH.main.but.statmaps,'Value',1);
  DATA.stim.curstimmodel = S.rgrs(1,:);
  showStatMaps;

  set(GH.main.but.threeDDisplay,'Value',1);
  start3DDisplay;

  % setup functional overlay
  set(GH.main.but.funcovly,'Value',1);
  setFunctionalOverlay;

  
function fillGlobalFieldsNeural(S,MCD) %,iChan)

  global DATA DSP

  % set DATA parameters explicitly for the combined data
  %
  voltc  = MCD.mvoltc{1}; %MCD.catmvoltc{1}; %MCD.catmvoltc;
  dgz    = MCD.catmdgz;
  imgdat = MCD.catimgdat;
  unimov = MCD.unimov;
  
  ncells = length(S);
  
%   for i=1:ncells
%      rgrs(i,:) = S(i).rgrsMION_resample'; %S(i).rgrsMION_resample';
%   end
% 
%   DATA.stim.curstattype  = 'rgr'; 
%   DATA.stim.stimmodel    = rgrs;
%   DATA.stim.curstimmodel = rgrs(iChan,:);
%   DATA.stim.curmodelname = 'REGRESSOR';
%   DSP.proc.corrmap_3d = S(iChan).mapR; %S(1).mapR;

  DSP.proc.params3d.nr    = size(voltc,4);
  DSP.proc.params3d.tr    = 2.4;
  [TR NR] = getTRandNR;
  DATA.ftimes       = [1:NR]*TR;  %[1:DSP.proc.params3d.nr]*DSP.proc.params3d.tr; %[1:NR]*TR;   
  DATA.valtp        = ones(1,length(DATA.ftimes));
  
  [lbtc,xs,ys]      = SU_convertVolsToLightboxes(voltc);
  DSP.lb.tc         = lbtc;        % only one long time course
  DSP.lb.anatvol    = SU_convertVolsToLightboxes(nanmean(voltc,4)); 
  DSP.lb.xs         = xs;
  DSP.lb.ys         = ys;

  DSP.lb.epivol     = nanmean(lbtc,3);
  DSP.lb.bkgvol     = nanmean(lbtc,3);
  DSP.proc.anatparams3d = DSP.proc.params3d;
  DSP.proc.anatvol_3d  = nanmean(voltc,4);
  DSP.proc.epivol_3d   = nanmean(voltc,4);
 
  DSP.proc.fmri_tc_3d= voltc; %%
  DATA.imagedat     = imgdat;
  DATA.dgz          = dgz;       
  
   
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MCD = loadMeanConcatData(textfileroot,skip,redo);
  
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
  
