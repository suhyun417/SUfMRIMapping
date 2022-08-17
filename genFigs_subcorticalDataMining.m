% genFigs_subcorticalDataMining.m
%
% 2022/06/16 SHP
% Playground to investigate face cells functional correlation with
% subcortical regions
%   - possible link to Daniel's resting state findings
%   - possible focus: known visual subcortical regions (pulvinar, sc, lgn),
%   basal forebrain, midbrain neuromodulatory centers, cerebellum,
%   claustrum (it's cortical but still)

clear all;

%% Settings
flagBiowulf = 0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/NIFVAULT/PROJECTS';
        directory.procdata = '/Volumes/NIFVAULT/PROCDATA';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/Volumes/NIFVAULT/LIBRARY';
        addpath(fullfile(directory.library, 'matlab_utils'));
    else % on virtual machine
        directory.projects = '/nifvault/projects';
        directory.procdata = '/nifvault/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/nifvault/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end


%% First up, are they driven by videos?
% split-half correlation across repeated viewings of videos
% excluding later sessions with MION b/c of the potential dropout in these
% regions

% let's load and compute correlation, make map and see
setSubjBOLD = {'Art', 'Ava'};

% for iSubj = 1:length(setSubjBOLD)
iSubj = 1;
nameSubjBOLD = setSubjBOLD{iSubj};
load(fullfile(directory.dataHome, ...
    sprintf('%s/_laterSessionExcluded/%s_movieTS_fMRI_indMov_splitHalf.mat', nameSubjBOLD, nameSubjBOLD)));

catTs_half = cell(2,1);
for iHalf = 1:2
    fmritc = [];
    for iM = 1:3
        curvoltc = voltcIndMov{iM, iHalf};
        avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
        if ~isempty(find(avgvoltc==0, 1))
            avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
        end
        pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
        fmritc = cat(4,fmritc,pcvoltc);
    end
    catTs_half{iHalf} = fmritc;
end

[nx, ny, nz, nt] = size(catTs_half{1});
nVox = nx*ny*nz;

catTs_half1_rs = reshape(catTs_half{1}, nVox, nt)';
catTs_half2_rs = reshape(catTs_half{2}, nVox, nt)';

matR_fMRIsplitHalf = NaN(nVox, 1);
for iVox = 1:nVox
    Rval = corr(catTs_half1_rs(:, iVox), catTs_half2_rs(:, iVox),...
        'rows','complete', 'type', 'Spearman');
    matR_fMRIsplitHalf(iVox, 1) = Rval;
end

matR_fMRIsplitHalf_3D = reshape(matR_fMRIsplitHalf, [nx, ny, nz]);

% % one-time save
% save(fullfile(directory.dataHome, ...
%     sprintf('%s/_laterSessionExcluded/%s_movieTS_fMRI_indMov_splitHalf.mat', nameSubjBOLD, nameSubjBOLD)), ...
%     'matR*', '-append')
% save(fullfile(directory.dataHome, ...
%     sprintf('%s/%s_movie123_splitHalfCorr_earlierSession.mat', nameSubjBOLD, nameSubjBOLD)), ...
%     'paramBOLD', 'voltcIndMov', 'matR*')

%% save a BRIK/HEAD registered to anatomy
dirDataBOLD = fullfile(directory.dataHome, nameSubjBOLD);
pname = [dirDataBOLD, '/tempSURF/'];
dirSPEC = [dirDataBOLD, '/Anatomy/_suma/'];

cd /nifvault/projects/parksh/_toolbox/BlockAna/
blockana;
S_neuralRegressor('Tor', nameSubjBOLD);
global STDPATH DSP DATA GH

for iMask = 1:2
    cd /nifvault/projects/parksh/_toolbox/BlockAna/
    
    switch iMask
        case 1 % no mask
            DSP.proc.fncvol_3d = matR_fMRIsplitHalf_3D; %.*brainMask_BlockAna3D;
            fname_tail = '';
        case 2 % brainmask
            load(fullfile(directory.dataHome, sprintf('/%s/%s_MaskArrays.mat', nameSubjBOLD, nameSubjBOLD)),...
                'brainMask_BlockAna3D');
            DSP.proc.fncvol_3d = matR_fMRIsplitHalf_3D.*brainMask_BlockAna3D;
            fname_tail = '_brainMask';
    end
    
    fname = sprintf('%s_Movie123_SplitHalfCorr%s+orig.BRIK', nameSubjBOLD, fname_tail);%
    
    vol = single(DSP.proc.fncvol_3d);
    
    % Do dumpFunctionalBrik
    Info = DATA.Info;
    Info.DATASET_RANK = DATA.Info.DATASET_RANK;
    Info.DATASET_DIMENSIONS = DATA.Info.DATASET_DIMENSIONS;
    Info.DELTA = DATA.Info.DELTA;
    Info.BYTEORDER_STRING = DATA.Info.BYTEORDER_STRING;
    Info.TYPESTRING = DATA.Info.TYPESTRING;
    Info.SCENE_DATA = DATA.Info.SCENE_DATA;
    Info.ORIGIN = DATA.Info.ORIGIN;
    Info.TAXIS_FLOATS = DATA.Info.TAXIS_FLOATS;
    Info.TAXIS_OFFSETS = DATA.Info.TAXIS_OFFSETS;
    Info.IDCODE_STRING = DATA.Info.IDCODE_STRING;
    Info.IDCODE_DATE = DATA.Info.IDCODE_DATE;
    Info.HISTORY_NOTE = DATA.Info.HISTORY_NOTE;
    
    Info.DATASET_DIMENSIONS = [size(vol) 0 0];
    Info.DATASET_RANK(2) = 1;
    Info.BRICK_TYPES = 3;  %3 = float
    Info.BRICK_FLOAT_FACS = 0;
    Info.BRICK_STATS = [min(vol(:)) max(vol(:))];
    Info.BRICK_LABS = DATA.Info.BRICK_LABS(1:2);
    Info.BRICK_KEYWORDS = DATA.Info.BRICK_KEYWORDS;
    Info.TAXIS_NUMS(1) = 1;
    Info.TypeName = 'float';
    Info.TypeBytes = 4;
    Info.ByteOrder = 'ieee-le';
    Info.Orientation = DATA.Info.Orientation;
    Info.FileFormat = DATA.Info.FileFormat;
    Info.RootName = DATA.Info.RootName;
    Info.Extension_1D = DATA.Info.Extension_1D;
    
    code0 = getAFNIOrientationCode('LR');
    code1 = getAFNIOrientationCode('AP');
    code2 = getAFNIOrientationCode('DV');
    
    Info.ORIENT_SPECIFIC = [code0 code1 code2];
    
    doti = strfind(fname,'+');
    prefix = fname(1:doti-1);
    Opt.Prefix = prefix;
    Opt.NoCheck = 0;
    
    fprintf('Writing %s as BRIK and HEAD files to %s\n',prefix,pname);
    curdir = cd(pname);
    
    %
    % Compensating for an AFNI inconsistency
    %
    if isfield(Info,'WARP_TYPE'),
        Info = rmfield(Info,'WARP_TYPE');
    end
    
    WriteBrik(vol,Info,Opt);  % needs to have a modified version to overwrite
    % an existing file with "w+"
    
    %%%% ADDED BY BER 4/5/13 for alignment purposes.
    % copy data to new file so original dimensions are saved.
    copyDataCMD = ['3dcopy  ' prefix ' ' prefix '_refit'];
    disp(copyDataCMD);
    system(copyDataCMD);
    
    % refit the x y z dimensioned based on strange AFNI conventions
    refitDataCMD = ['3drefit -xdel 1.5 -ydel 1.5 -zdel 1.5 ' prefix '_refit+orig'];
    disp(refitDataCMD);
    system(refitDataCMD);
    
    monk=DSP.afni_prefix(1:3);
    % now that data is refit, align the center with the surface anatomy center
    alignDataCMD = ['@Align_Centers -base ' STDPATH.afni_uw '/Anatomy/' monk '_t1_refit_shft+orig. -dset ' prefix '_refit+orig.'];
    disp(alignDataCMD);
    system(alignDataCMD);
    
    % finally shift the origin slightly because of strange ANFI stuff
    shiftDataCMD = ['3drefit -dyorigin -0.5 -dzorigin 0.5 ' prefix '_refit_shft+orig.'];
    disp(shiftDataCMD);
    system(shiftDataCMD);
    %%%%
    
    % Copy relevant files to where the surface files are
    [s,mess,messID]=movefile([pname, '*_refit_shft+orig*'],  dirSPEC);
    if ~s
        fprintf(1,' %s \n', mess);
    end
end
    
    %% save split-half averaged data
    % setSubjBOLD = {'Art', 'Ava'};
    %
    % for iSubj = 1:length(setSubjBOLD)
    %     nameSubjBOLD = setSubjBOLD{iSubj};
    %     saveFileName = fullfile(directory.dataHome, ...
    %         sprintf('%s/_laterSessionExcluded/%s_movieTS_fMRI_indMov_splitHalf.mat', nameSubjBOLD, nameSubjBOLD));
    %
    %     % load(fullfile(directory.dataHome, sprintf('%s/_laterSessionExcluded/%s_movieTS_fMRI_indMov.mat', nameSubjBOLD, nameSubjBOLD)));
    %
    %     % Get information from the filelist in .txt
    %     switch lower(nameSubjBOLD)
    %         case 'art'
    %             sessionFileList = 'FL_e66_allfiles2.txt';
    %         case 'ava'
    %             sessionFileList = 'FL_ava_allfiles2.txt';
    %     end
    %
    %     fprintf(1, '\n Filelist in *.txt is "%s" \n', sessionFileList);
    %
    %     skip=7; % number of TRs to skip
    %
    %     fprintf(1, '\n Skipping first %d TRs of each movie by replacing them with NaNs \n', skip);
    %
    %     % cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
    %
    %     % S_neuralRegressor(nameSubjNeural, nameSubjBOLD);
    %     % global STDPATH DSP DATA GH
    %
    %     addpath('/nifvault/projects/parksh/_toolbox/BlockAna/')
    %     blockana;
    %     global GH STDPATH DSP
    %
    %     SI = SU_createSessInfo(sessionFileList,[],skip);
    %     filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);
    %
    %     MAX_TR = SI.max_tr;
    %     skip=SI.skiptr;
    %
    %     % unimov = sort(unique(SI.movID));
    %     setMovie = [1 2 3];
    %
    %     % nunimov = length(unimov);
    %     catimgdat = [];
    %     % rgr = zeros(1,nunimov*MAX_TR);
    %     % procmovs = [];
    %     infoSession = struct([]);
    %     for iMovie = 1:length(setMovie)
    %         movindxs = find(SI.movID == setMovie(iMovie));
    %         nmovindxs = length(movindxs);
    %
    %         % collect all the data files for this movie
    %         %
    %         filelist = SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs),...
    %             SI.scanID(movindxs));
    %
    %         totfiles = length(filelist);
    %
    %         %% get the information on acquisition date and exclude later sessions
    %         clear setDates indValidSession
    %         for f=1:totfiles
    %             s_sub = filelist{f}.subj;
    %             s_ses = filelist{f}.sess;
    %             s_sc  = filelist{f}.scan;
    %
    %             cursession = sprintf('%s.%s',s_sub,s_ses);
    %             sessindx   = strmatch(cursession,DSP.filedetails.sess);
    %             scanindx   = find(DSP.filedetails.scan == s_sc);
    %             xlsindx    = intersect(sessindx,scanindx);
    %             setDates{f, 1} = DSP.filedetails.date{xlsindx};
    %         end
    %         critDate = '01/01/12';
    %         critDateNum = datenum(critDate, 'mm/dd/yy');
    %         indValidSession = find(datenum(string(setDates)) < critDateNum);
    %
    %         infoSession(iMovie).nameSubj = nameSubjBOLD;
    %         infoSession(iMovie).movieID = setMovie(iMovie);
    %         infoSession(iMovie).setSessionDates = setDates;
    %         infoSession(iMovie).critDate = critDate;
    %         infoSession(iMovie).indValidSession = indValidSession;
    %         infoSession(iMovie).setSessionDates_valid = setDates(indValidSession);
    %
    %
    %         %%
    %         filelist_valid = filelist(indValidSession); %SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs), SI.scanID(movindxs));
    %
    %         totfiles_valid = length(filelist_valid);
    %
    %         for f=1:totfiles_valid
    %             s_sub = filelist_valid{f}.subj;
    %             s_ses = filelist_valid{f}.sess;
    %             s_sc  = filelist_valid{f}.scan;
    %
    %             [fmri_tc,dgz,notes] = SU_loadFile(s_sub,s_ses,s_sc);
    %
    %             if f==1
    %                 [a,b,c,t] = size(fmri_tc);
    %                 allvoltc = zeros(a,b,c,MAX_TR,totfiles);
    %             end
    %
    %             scanlen = size(fmri_tc,4);
    %             if scanlen == 250
    %                 mov_start_tr = 63;% offset when movie starts in scan
    %             elseif scanlen == 125
    %                 mov_start_tr = 1; % short movie
    %             else
    %                 fprintf('WARNING: bad scan length %d\n',scanlen);
    %             end
    %
    %             val_tr         = [mov_start_tr:(mov_start_tr+MAX_TR-1)];
    %             [clp_fmri_tc,clp_mdgz] = SU_clipMovDat(fmri_tc,dgz,val_tr);
    %             allvoltc(:,:,:,:,f) = clp_fmri_tc;
    %             alldgz(f) = dgz;
    %         end
    %
    %         % take half and average separately
    %         clear orderFiles
    %         orderFiles = randperm(totfiles_valid);
    %         if rem(totfiles_valid, 2) % odd number
    %             orderFiles = orderFiles(1:end-1);
    %         end
    %         orderFiles = reshape(orderFiles, 2, length(orderFiles)/2);
    %
    %         infoSession(iMovie).filelist_valid = filelist_valid;
    %         infoSession(iMovie).splithalf_indFiles = orderFiles;
    %
    %         for iHalf = 1:2
    %
    %             [mvoltc,mdgz]=SU_computeMeanTC(allvoltc(:,:,:,:,orderFiles(iHalf, :)), alldgz(orderFiles(iHalf, :)));
    %
    %             if skip ~= 0
    %                 vox_mns=nanmean(mvoltc,4);
    %                 for i=1:skip
    %                     % sets the skipped TRs to be the mean of the voxel 10/30/12 BER
    %                     %
    %                     mvoltc(:,:,:,i)=vox_mns(:,:,:);
    %                 end
    %             end
    %
    %             fprintf('Starting at vol %d. Detrending mov %d\n',mov_start_tr, setMovie(iMovie));
    %             dt_mvoltc = SU_detrendMovDat(mvoltc); %mvoltc;
    %
    %             if skip ~=0
    %                 % sets the skipped TRs to be NAN 10/30/12 BER
    %                 %
    %                 dt_mvoltc(:,:,:,1:skip)=NaN;
    %             end
    %
    %
    %             %     allmvoltc(:,:,:,:,iMovie) = dt_mvoltc;
    %             %     allmdgz(iMovie) = clp_mdgz;
    %
    %             voltcIndMov{iMovie, iHalf} = dt_mvoltc;
    %
    %             %     catimgdat = cat(2, catimgdat, SU_getImageDat(unimov(iMovie)));
    %             %     end
    %
    %         end
    %
    %         % parameters
    %         paramBOLD.infoSession = infoSession;
    %         paramBOLD.max_tr = MAX_TR;
    %         paramBOLD.Fs = SI.Fs;
    %         paramBOLD.TR = SI.TR;
    %         paramBOLD.skiptr = SI.skiptr;
    %         paramBOLD.unimov = setMovie; %unimov;
    %
    %         fprintf(1, ':::: %s: Movie %d: # of Trials %d for half :::: \n', ...
    %             nameSubjBOLD, setMovie(iMovie), length(orderFiles))
    %
    %
    %     end
    %
    %     % save the data file
    %     save(saveFileName, 'voltcIndMov', 'paramBOLD'); %, 'catData');
    %     fprintf(1, 'Time series are saved in %s\n', saveFileName)
    %
    % end
    %
