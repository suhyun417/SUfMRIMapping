function [] = saveMovieDataBOLD_indMov() %, saveFileName_matSDF, saveFileName) %[dataBOLD] = genMovieDataBOLD_indMov(nameSubjBOLD) %,
% 
% Save fMRI data for individual movies for a particular subject under the individual data folder
% Modified from "genMovieDataBOLD_indMov.m" and "savemovietseries.m"
% 3/20/15, SHP

% addpath ('/projects/parksh')

% Subject name and directory
fprintf(1, '\nSave fMRI time series for individual movies separately\n');
nameSubjBOLD = input('\nName of subject? (e.g. Art, Ava, Sig):', 's');
if sum(strcmpi(nameSubjBOLD, {'artemis', 'art', 'ar', 'e66'}))
    nameSubjBOLD = 'Art';
elseif sum(strcmpi(nameSubjBOLD, {'avalanche', 'ava', 'av'}))
    nameSubjBOLD = 'Ava';
elseif sum(strcmpi(nameSubjBOLD, {'Ziggy', 'Sieglinde', 'Si', 'Sieggie', 'Sig', 'AZ26', 'AZ2'}))
    nameSubjBOLD = 'Sig';
end

dirData = '/procdata/parksh/macaque';
dirDataBOLD = [dirData, nameSubjBOLD, '/'];
saveFileName = fullfile(dirDataBOLD, sprintf('%s_movieTS_fMRI_indMov.mat', nameSubjBOLD));

fprintf(1, '\n Save fMRI time series for subject "%s"', nameSubjBOLD);
fprintf(1, '\n Data will be saved as "%s" \n ', saveFileName);
fprintf(1, '\nPress enter to proceed \n'); input('')

% Get information from the filelist in .txt
switch lower(nameSubjBOLD)
    case 'art'
        sessionFileList = 'FL_e66_allfiles2.txt';
    case 'ava'
        sessionFileList = 'FL_ava_allfiles2.txt';
    case 'sig'
        sessionFileList = 'FL_sig_allfiles.txt';
end

fprintf(1, '\n Filelist in *.txt is "%s" \n', sessionFileList);

skip=7; % number of TRs to skip

fprintf(1, '\n Skipping first %d TRs of each movie by replacing them with NaNs \n', skip);

% cd /projects/parksh/NeuralBOLD/analysis/BlockAna/

% S_neuralRegressor(nameSubjNeural, nameSubjBOLD);
% global STDPATH DSP DATA GH

addpath('/projects/parksh/_toolbox/BlockAna/')
blockana;

SI = SU_createSessInfo(sessionFileList,[],skip); 
filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);

% % Init by loading one file
% s_sub = filelist{1}.subj;
% s_ses = filelist{1}.sess;
% s_sc  = filelist{1}.scan;
% [fmri_tc_3d,dgz] = SU_loadFile(s_sub,s_ses,s_sc);
% 
% [TR,NR] = getTRandNR();

% if ~isempty(dir(dirDataBOLD)) % & ~redo
%     load(sessdatafile);  % load already created file
% else
    % Create and save the mean tseries for each movie
    %
%     dataBOLD = buildTimeSeries(SI);     % create the big S structure
             % save the file
% end


% MCD = loadMeanConcatData(sessionFileList,skip,redo); 
% function MCD = loadMeanConcatData(textfileroot,skip,redo)
%   

% function MCD = buildTimeSeries(SI)

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
                        
        fmri_tc_3d = DSP.proc.fmri_tc_3d;
        dgz = [];
        if isfield(DATA,'dgz')
            dgz     = DATA.dgz;
        end
        
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
    dt_mvoltc = SU_detrendMovDat(mvoltc); %mvoltc;
    
    if skip ~=0
        % sets the skipped TRs to be NAN 10/30/12 BER
        %
        dt_mvoltc(:,:,:,1:skip)=NaN;
    end
    
    
    allmvoltc(:,:,:,:,u) = dt_mvoltc;
    allmdgz(u) = clp_mdgz;
    
    voltcIndMov{u} = dt_mvoltc;
    
    catimgdat = cat(2,catimgdat,SU_getImageDat(unimov(u)));
end

[catmvoltc,catmdgz]  =...
    SU_concatVolumesAndDGZ(allmvoltc,allmdgz);
catData.catimgdat = catimgdat;
catData.catmdgz = catmdgz;


% parameters
paramBOLD.max_tr = MAX_TR;
paramBOLD.Fs = SI.Fs;
paramBOLD.TR = SI.TR;
paramBOLD.skiptr = SI.skiptr;
paramBOLD.unimov = unimov;

% [catmvoltc,catmdgz]  =...
%     SU_concatVolumesAndDGZ(allmvoltc,allmdgz);

% Store Scan Info in S % Added by BER 02/12/13
%
% MCD.max_tr=MAX_TR;
% MCD.Fs=SI.Fs;
% MCD.TR=SI.TR;
% MCD.skiptr=SI.skiptr;
% 
% MCD.catmvoltc = catmvoltc;
% MCD.catmdgz   = catmdgz;
% MCD.catimgdat = catimgdat;
% MCD.unimov    = unimov;

% save the data file
save(saveFileName, 'voltcIndMov', 'paramBOLD', 'catData');  
fprintf(1, 'Time series are saved in %s\n', saveFileName)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% hier 


% % Get the movie list from the fMRI data
% fprintf(1, '\Loading BOLD time series...\n');
% 
% d_b = dir(fullfile(dirDataBOLD, '*BOLD.mat')); %W = what(dirDataBOLD);
% load(fullfile(dirDataBOLD, d_b(1).name)); % for now it's only one file
% setMovIDs = MCD.unimov; % Movie IDs of the current BOLD
% 
% 
% %% First, split fMRI responses into each individual movie
% 
% fprintf(1, 'Number of movies in original time series: %d\n', length(setMovIDs))
% 
% dataBOLD = MCD; % Get necessary params for further DSP in BlockAna: for now, just all the params
% dataBOLD = rmfield(dataBOLD, 'catmvoltc'); %{'catmvoltc', 'catmdgz', 'catimgdat'}); % remove the previous variable of concatenated time series
% 
% lengthIndMovie_TR = dataBOLD.max_tr;
% % lengthIndMovie_sec = lengthIndMovie_TR*dataBOLD.TR; % 300;
% 
% % Get the necessary part of time courses, corresponding to our (common) movie set
% % if ~issorted(setMovIDs) % just to be sure
% %     setMovIDs = sort(setMovIDs);
% % end
% 
% % countTR = 0; tc=[];
% % Split concatenated time series to individual cells for individual movies
% for iMov = 1:length(setMovIDs)
%     fprintf(1, 'Movie ID: %d (%d/%d), Location in original BOLD: %d\n', setMovIDs(iMov), iMov, length(setMovIDs), setMovIDs(iMov));
%     %     tc(:,:,:,countTR+1:countTR+lengthIndMovie_TR) = MCD.catmvoltc(:,:,:,lengthIndMovie_TR*(locMovBOLD(iMov)-1)+1:lengthIndMovie_TR*locMovBOLD(iMov)); % should be optimized for different set of movies
%     dataBOLD.mvoltc{iMov} = MCD.catmvoltc(:,:,:,lengthIndMovie_TR*(setMovIDs(iMov)-1)+1:lengthIndMovie_TR*setMovIDs(iMov)); % should be optimized for different set of movies
%     %     countTR = size(tc, 4);
% end
% % dataBOLD.catmvoltc = tc;
% % fprintf(1, 'Time course for %d movies has been saved in dataBOLD.tc\n', length(setMovIDs));
% % fprintf(1, 'The number of TRs of dataBOLD.catmvoltc should be %d and is %d\n', ...
% %     lengthIndMovie_TR*length(setMovIDs), size(dataBOLD.catmvoltc, 4));
% 
% dataBOLD.unimov = setMovIDs;
% 
% clear MCD %tc %SigData
% 
% % save File
% save(fullfile(dirDataBOLD, saveFileName), 'dataBOLD')


% % Select common movie IDs between single units and BOLD (or define the movies)
% setMovIDs = input('Movie IDs (e.g. [1 2 3])? (default: common movies in single units and fMRI)');
% if isempty(setMovIDs)
%     setMovIDs = intersect(listMovBOLD, str2num(char(listMovSU))');
% end
% % setMovIDs = [13 14 15]; %intersect(listMovBOLD, str2num(char(listMovSU))');
% [~, locMovBOLD] = intersect(listMovBOLD, setMovIDs);
% % [setMovIDs, locMovBOLD] = intersect(listMovBOLD, str2num(char(listMovSU))');
%
% % flagConcat = input('Do you want to concatenate timeseries across movies? (1:yes, 0:no):');
%
% % Get the list of single units & their movie
% listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
% listSUchannelID = unique(listSU_all);
%
% indDataMov=[];
% for iCh=1:length(listSUchannelID) % go throucgh channel-by-channel
%      tempListMov{iCh} = regexp([listMatSUFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
%      % Get indices for common movies across cells
%      indDataMov(iCh, :) = ismember(setMovIDs, str2num(char(tempListMov{iCh}))'); % data presence matrix (Channels x MovieIDs)
% end
%
% % % For now, let's do it with cells that have data for every movies
% % setCellIDs = listSUchannelID(sum(indMov,2)==length(setMovIDs));
% setCellIDs = listSUchannelID;
%
% % setCellIDs = {'006a'; '009a'; '010a'; '012a'; '013a'}; % Rhombus %{'003a';'005a';'007a';'009a';'013b';'014a'}; % Sig: '013a' has very low spikes
% % setMovIDs  = [7 8 9 10 11 12]; % Rhom %[7 8 9]; % Sig % [1 2 3 7 8 9];





