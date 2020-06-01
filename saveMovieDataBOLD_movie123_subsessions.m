% saveMovieDataBOLD_movie123_subsessions.m
% Script;
%
% 2020/05/11 SHP: Exclude some later sessions that show some dropout due to
% MION accumulation. Previous file including all the sessions are named now
% "Art_movieTS_fMRI_indMov_all.mat"
%
% Save fMRI data for individual movies for a particular subject under the individual data folder
% Modified from "genMovieDataBOLD_indMov.m" and "savemovietseries.m"
% 3/20/15, SHP

% addpath ('/projects/parksh')

setNameSubj = {'Art', 'Ava'};
for iSubj = 1:2
    nameSubjBOLD = setNameSubj{iSubj};
    
    dirData = '/procdata/parksh/_macaque';
    dirDataBOLD = fullfile(dirData, nameSubjBOLD);
    saveFileName = fullfile(dirDataBOLD, sprintf('%s_movieTS_fMRI_indMov.mat', nameSubjBOLD));
    
    fprintf(1, '\n Save fMRI time series for subject "%s"', nameSubjBOLD);
    fprintf(1, '\n Data will be saved as "%s" \n ', saveFileName);
    % fprintf(1, '\nPress enter to proceed \n'); input('')
    
    % Get information from the filelist in .txt
    switch lower(nameSubjBOLD)
        case 'art'
            sessionFileList = 'FL_e66_allfiles2.txt';
        case 'ava'
            sessionFileList = 'FL_ava_allfiles2.txt';
    end
    
    fprintf(1, '\n Filelist in *.txt is "%s" \n', sessionFileList);
    
    skip=7; % number of TRs to skip
    
    fprintf(1, '\n Skipping first %d TRs of each movie by replacing them with NaNs \n', skip);
    
    % cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
    
    % S_neuralRegressor(nameSubjNeural, nameSubjBOLD);
    % global STDPATH DSP DATA GH
    
    addpath('/projects/parksh/_toolbox/BlockAna/')
    blockana;
    global GH STDPATH DSP
    
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
    
    % unimov = sort(unique(SI.movID));
    setMovie = [1 2 3];
    
    % nunimov = length(unimov);
    catimgdat = [];
    % rgr = zeros(1,nunimov*MAX_TR);
    % procmovs = [];
    infoSession = struct([]);
    for iMovie = 1:length(setMovie)
        movindxs = find(SI.movID == setMovie(iMovie));
        nmovindxs = length(movindxs);
        
        % collect all the data files for this movie
        %
        filelist = SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs),...
            SI.scanID(movindxs));
        
        totfiles = length(filelist);
        
        %% get the information on acquisition date and exclude later sessions
        clear setDates indValidSession
        for f=1:totfiles
            s_sub = filelist{f}.subj;
            s_ses = filelist{f}.sess;
            s_sc  = filelist{f}.scan;
            
            cursession = sprintf('%s.%s',s_sub,s_ses);
            sessindx   = strmatch(cursession,DSP.filedetails.sess);
            scanindx   = find(DSP.filedetails.scan == s_sc);
            xlsindx    = intersect(sessindx,scanindx);
            setDates{f, 1} = DSP.filedetails.date{xlsindx};
        end
        critDate = '01/01/12';
        critDateNum = datenum(critDate, 'mm/dd/yy');
        indValidSession = find(datenum(string(setDates)) < critDateNum);
        
        infoSession(iMovie).nameSubj = nameSubjBOLD;
        infoSession(iMovie).movieID = setMovie(iMovie);
        infoSession(iMovie).setSessionDates = setDates;
        infoSession(iMovie).critDate = critDate;
        infoSession(iMovie).indValidSession = indValidSession;
        infoSession(iMovie).setSessionDates_valid = setDates(indValidSession);
        
        
        %%
        filelist_valid = filelist(indValidSession); %SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs), SI.scanID(movindxs));
        
        totfiles_valid = length(filelist_valid);
        
        for f=1:totfiles_valid
            s_sub = filelist_valid{f}.subj;
            s_ses = filelist_valid{f}.sess;
            s_sc  = filelist_valid{f}.scan;
            
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
        
        fprintf('Starting at vol %d. Detrending mov %d\n',mov_start_tr, setMovie(iMovie));
        dt_mvoltc = SU_detrendMovDat(mvoltc); %mvoltc;
        
        if skip ~=0
            % sets the skipped TRs to be NAN 10/30/12 BER
            %
            dt_mvoltc(:,:,:,1:skip)=NaN;
        end
        
        
        %     allmvoltc(:,:,:,:,iMovie) = dt_mvoltc;
        %     allmdgz(iMovie) = clp_mdgz;
        
        voltcIndMov{iMovie} = dt_mvoltc;
        
        %     catimgdat = cat(2, catimgdat, SU_getImageDat(unimov(iMovie)));
    end
    
    % [catmvoltc,catmdgz]  =...
    %     SU_concatVolumesAndDGZ(allmvoltc,allmdgz);
    % catData.catimgdat = catimgdat;
    % catData.catmdgz = catmdgz;
    
    
    % parameters
    paramBOLD.infoSession = infoSession;
    paramBOLD.max_tr = MAX_TR;
    paramBOLD.Fs = SI.Fs;
    paramBOLD.TR = SI.TR;
    paramBOLD.skiptr = SI.skiptr;
    paramBOLD.unimov = setMovie; %unimov;
    
    
    % save the data file
    save(saveFileName, 'voltcIndMov', 'paramBOLD'); %, 'catData');
    fprintf(1, 'Time series are saved in %s\n', saveFileName)
    
end


