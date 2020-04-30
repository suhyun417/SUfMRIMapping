function S_correlate_ALLR_sig()
  
  global S
  %
  % examination of regressors combined across days
  %
  
  global STDPATH DSP DATA GH
  
  
  % first start,open the XL file, and reset the windows
  %
  blockpaths_ber;
  STDPATH.xlsfile = 'block_timecourse_files_ber.xls';
  loadXLSFile;

  resetBlockanaWindows;


  % save as data file in session directory 
  %
  saveflag = 0;  loadS = 1; %9/12/13
  S_filename = 'sig_S_rgr_corr_lots';
  
  redo = 0;
  skip=7; % number of TRs to skip  %% ADDED 10/04/12 BER

  textfileroot = 'FL_sig_allfiles';
%   textfileroot = 'e66_2file_allmovies';

  MCD  = loadMeanConcatData(textfileroot,skip,redo); % load or

  % 
  % here read in some pre-computed regressors and read some from the xls
  % file (all in raw form)
  %
  RAWRGR = getAllRawRegressors(MCD.unimov,MCD.TR,MCD.Fs,MCD.max_tr);

  %
  % Here we replace the initial list of regressors 
  % with some that we tailor (e.g. combining, compressing, etc)
  %
%   ALLRGR = SU_addSpecificRegressors_Motion(RAWRGR,MCD.TR); 
%   ALLRGR = SU_addSpecificRegressors_Motion2(RAWRGR,MCD.TR); 
%   ALLRGR = SU_addSpecificRegressors_Faces(RAWRGR,MCD.TR); 
%   ALLRGR = SU_addSpecificRegressors_LowVision(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_Bodies(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_Behav(RAWRGR,MCD.TR);
  ALLRGR = SU_addSpecificRegressors_Compare(RAWRGR,MCD.TR);
%   ALLRGR = SU_addSpecificRegressors_ShortCompare(RAWRGR,MCD.TR);

  %
  % compute the correlation maps (single variable regression)
  %
  %
  %  ORDER OF FEATURES
    %  1   =  Speed(log_MION)
    %  2   =  Mot_Diverg_Rect(log_MION)
    %  3   =  Mot_Diverg_STD(log_MION)
    %  4   =  Mot_Diverg_Beta(log_MION)
    %  5   =  Mot_Diverg_Gamma(MION_sqrt_sm2)
    %  6   =  Scene_Cuts(log_MION_sm2)
    %  7   =  Faces_Full(sqrt_MION)
    %  8   =  Faces_SV(sqrt_MION)
    %  9   =  1_Face_Full(sqrt_MION)
    %  10  =  1_Face_SV(sqrt_MION)
    %  11  =  Heads(sqrt_MION)
    %  12  =  Faces_all(sqrt_MION)
    %  13  =  1_Face_comb(sqrt_MION)
    %  14  =  Luminance(log_MION)
    %  15  =  Contrast_STD(log_MION)
    %  16  =  Contrast_Beta(log_MION)
    %  17  =  Contrast_Gamma(log_MION)
    %  18  =  SF_low(log_MION)
    %  19  =  SF_high(log_MION)
    %  20  =  SF_ratio(log_MION)
    %  21  =  Butts(MION_sqrt)
    %  22  =  Bodies(sqrt_MION)
    %  23  =  Extremities(sqrt_MION)
    %  24  =  Hands(sqrt_MION)
    %  25  =  Animals(MION_sqrt)
    %  26  =  Macaque(MION_sqrt)
    %  27  =  Rhesus(MION)
    %  28  =  Human(MION)
    %  29  =  Mot_Local_STD(MION)
    %  30  =  Mot_Bartel_Local(MION)
    %  31  =  Mot_Bartel_Global(MION)
    %  32  =  Mot_Bartel_Resid(MION)
    %  33  =  Mot_Bartel_Total(MION)
    
  if isempty(S) && loadS == 0

    rgrlist=[1 29 13 12 23 25 16 17 14 ]; % [1:Global_Speed 2:Local_Motion 3:1_Face_comb 4:Faces_all ....
                                   %  5:Extremities 6:Animals 7:Contrast_Beta 8:Contrast_Gamma 9:Luminance]
      
    ALLRGR.rgrs=ALLRGR.rgrs(:,rgrlist);
    ALLRGR.names=ALLRGR.names(rgrlist);
    [RMAP PMAP] = computeIndCorrMaps(MCD.catmvoltc,ALLRGR.rgrs);
  
    S.rgrs     = ALLRGR.rgrs;
    S.names    = ALLRGR.names;
    S.pmaps = PMAP;  % Need to change this (!)
    S.corrmaps  = (RMAP.*-1);  % no longer rsq map now just correlations but inverted because of MION
  elseif loadS == 1
    cd(STDPATH.sessdatdir)
    eval(['load(''' S_filename ''');'])
    cd(STDPATH.blockana_home)
  else
    disp('S preloaded in memory so skipping regression');
  end
  
  DSP.proc.fracvarmap_3d = squeeze(S.corrmaps(:,:,:,1));
      
  fillGlobalFields(MCD,S);
  showVolumeTimeCourse;

  % setup advanced functional
  set(GH.main.but.advanced,'Value',1);
  doAdvancedFunctional;
  
  set(GH.main.but.statmaps,'Value',1);
  DATA.stim.curstimmodel = S.rgrs(:,1);
  showStatMaps;

  set(GH.main.but.threeDDisplay,'Value',1);
  start3DDisplay;

  % setup functional overlay
  set(GH.main.but.funcovly,'Value',1);
  setFunctionalOverlay;
      
  if saveflag
      cd(STDPATH.sessdatdir)
      eval(['save ' S_filename ' S'])
      cd(STDPATH.blockana_home)
  end
  

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
  

function RAWRGR = getAllRawRegressors(unimov,TR,Fs,max_tr)
  codingpath = '/einstein0/USRlab/data/russbe/CodingXLS';
  lowlevpath = '/einstein0/USRlab/data/russbe/LowlevelMovRegressors';
  
  codingfile = sprintf('CodingSheetmov1_inv');
  xlsfullpath = sprintf('%s/%s.xls',codingpath,codingfile);
  [num,str] = xlsread(xlsfullpath);
  varnames = str(1,:);
 
  lowlevfile = sprintf('Movie1_rgr');
  lowlevfullpath = sprintf('%s/%s.mat',lowlevpath,lowlevfile);
  load(lowlevfullpath);  % loads RGR variable
  
  varnames = [varnames RGR.features];
    
  if nargin < 2  %% ADDED skipTR variable 10/04/12 BER
      skipTRs=0; % don't skip any TRs if no skip is sent
  end
  
  cregressors = [];
  rgrs = [];
  for m=1:length(unimov)
	codingfile   = sprintf('CodingSheetmov%d_inv',unimov(m));
	lowlevfile = sprintf('Movie%d_rgr',unimov(m));
	
	xlsfullpath = sprintf('%s/%s.xls',codingpath,codingfile);
	[num,str] = xlsread(xlsfullpath);
	lowlevfullpath = sprintf('%s/%s.mat',lowlevpath,lowlevfile);
	load(lowlevfullpath);  % loads RGR variable
	
	
	nvars = length(str(1,:));
	for i=1:nvars
	  tmp   = num(:,i);
	  rtmp  = resample(tmp,100*Fs,100*TR);
	  rtmp(max_tr+1:end) = [];
	  rgrs(:,i) = rtmp;
	end
	
	nrgr = size(RGR.regressors,1);
	for i=1:nrgr
	  rgrs(:,nvars+i)=RGR.regressors(i,:);
	end
	
% 	tvals = [(-20*TR):TR:(20*TR)];
% 	% This is just made-up
% 	MION_k = gampdf(tvals,TR,2*TR);
% 	MION_k = MION_k./sum(MION_k);
% 	crgrs = doConv(rgrs,MION_k);
% 	cregressors = [cregressors crgrs];  
    cregressors = [cregressors rgrs'];
  end
  
  RAWRGR.names = varnames;
  RAWRGR.rgrs  = cregressors';
  
  
function [Rval_maps Pval_maps] = computeIndCorrMaps(voltc,rgrs)

  [nx,ny,nz,nt] = size(voltc);
  nrgr = size(rgrs,2);
  voltc_lin = reshape(voltc,[nx*ny*nz nt])';  % make 2-D
  
  nvox = size(voltc_lin,2);
  Rvals_lin = zeros(nrgr,size(voltc_lin,2));
  Pvals_lin = zeros(nrgr,size(voltc_lin,2));
  
  wb = waitbar(0, 'Computing Regression');

  for r=1:nrgr
	waitbar(r/nrgr,wb)
	[tmpRvals, tmpPvals] = corr(voltc_lin,real(rgrs(:,r)),'rows','complete');
	Pvals_lin(r,:)=tmpPvals;
	Rvals_lin(r,:) = tmpRvals;
  end
  close(wb)

  Rval_maps = reshape(Rvals_lin',[nx ny nz nrgr]);
  Pval_maps = reshape(Pvals_lin',[nx ny nz nrgr]);
 
  
  
function fillGlobalFields(MCD,S)

  global DATA DSP

  % set DATA parameters explicitly for the combined data
  %
  voltc  = MCD.catmvoltc;
  dgz    = MCD.catmdgz;
  imgdat = MCD.catimgdat;
  unimov = MCD.unimov;
  
  DSP.proc.params3d.nr    = size(voltc,4);
  [TR NR] = getTRandNR;
  DATA.ftimes       = [1:NR]*TR;   
  DATA.valtp        = ones(1,length(DATA.ftimes));
  
  [lbtc,xs,ys]      = SU_convertVolsToLightboxes(voltc);
  DSP.lb.tc         = lbtc;        % only one long time course
  DSP.lb.anatvol = SU_convertVolsToLightboxes(nanmean(DSP.proc.fmri_tc_3d,4));
  DSP.lb.xs         = xs;
  DSP.lb.ys         = ys;

  DSP.lb.epivol     = nanmean(lbtc,3);
  DSP.lb.bkgvol     = nanmean(lbtc,3);
  DSP.proc.anatparams3d = DSP.proc.params3d;
  DSP.proc.anatvol_3d  = nanmean(voltc,4);
  DSP.proc.epivol_3d   = nanmean(voltc,4);
 
  DSP.proc.fmri_tc_3d= voltc;
  DATA.imagedat     = imgdat;
  DATA.dgz          = dgz;       
  
  DATA.stim.curstattype  = 'rgr'; 
  DATA.stim.stimmodel    = S.rgrs;
  DATA.stim.curstimmodel = S.rgrs(:,1);
  DATA.stim.curmodelname = 'REGRESSOR';

  showVolumeTimeCourse;
  setLightboxOrientation('hor');
  updateLightbox;
   
  
% function startFracVarAnalysis(MCD,S)
% 
%   global GH
%   
%   % load into blockana for viewing, with widget to switch between regressors
%   %
%   %
%   nrgr = size(S.rgrs,2);
%   
%   INIT = 1;
% 
%   figure(3020);clf;
%   set(gcf,'Units','Normalized','Position',[0.3 0.6 0.3 0.03])
%   GH.SU.rgrsel = uicontrol('Style','popup','Units','normalized',...
% 						   'Position',[0.05 .2 0.9 0.7],...
% 						   'String', S.names,...
% 						   'Callback',sprintf('SU_loadFracVarMap'),...
% 						   'UserData',S);
%     
%   fillGlobalFields(MCD,S);
%   showVolumeTimeCourse;
%   
%   % setup advanced functional
%   set(GH.main.but.advanced,'Value',1);
%   doAdvancedFunctional;
%   
%   % setup tcexplore
%   %set(GH.but.tcexplore,'Value',1);
%   %tcExplore;
%   
%   SU_loadFracVarMap(squeeze(S.betamaps(:,:,:,INIT)),...
% 					squeeze(S.rsqmaps(:,:,:,INIT)),S.rgrs(:,INIT));
%   
%   % setup 3d Display
%   set(GH.main.but.threeDDisplay,'Value',1);
%   start3DDisplay;
%   
%   % setup functional overlay
%   set(GH.main.but.funcovly,'Value',1);
%   setFunctionalOverlay;
%   