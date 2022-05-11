% saveMovieData_SU_EyeSignal_ForMarta.m
% 2022/04/14 SHP
% modified from "saveMovieData_SU_EyeSignal_ForReza.m"
%   - in addition to the movie data + eye movement signal from some
%   sessions, putting together the data from fingerprinting sessions


% tankdata.streams.Anlg.data : x , y, photodiode, reward, pupil diameter
%                                                   (from channel 1 to 5)
% tankdata.streams.Anlg.fs : analog signal sampling rate
% tankdata.scalars.Evnt.data : event code (100: the start of the movie for movie sesions, 101 - 1101: image ID for fingerprinting)
% tankdata.scalars.Evnt.ts: timestamp of event code. in seconds.
% celldata.ts: timestamp of spikes. in milliseconds.
% celldata.ts and tankdata.streams.Anlg, tankdata.scalars.Evnt.ts are
% aligned to the start of the recording

% cd /procdata/parksh/_macaque/_Marta

clear all;

% load transformation matrix
load /nifvault/projects/koyanok/Codes/BehaviorControl/MonkeyLogic/Volt2DegRig1.mat;

flagSaveMovie = 1; %0; % movie
flagSaveFP = 0; %1; % fingerprinting

setSession = {'Mochi181023', 'Mochi181024'}; %{'Mochi190313', 'Mochi190314'}; %{'Mochi181023', 'Mochi181024'}; %{'Wasabi190304', 'Wasabi190305', 'Wasabi190306'}; %{'Davida180723', 'Davida180724', 'Davida180725'}; %{'Spice180124', 'Spice180126'}; %

for iSession = 1:length(setSession)
    nameSession = setSession{iSession};
    
    if contains(nameSession, 'Spice')
        d_session = dir(sprintf('/nifvault/procdata/koyanok/physiology/cells/%s', nameSession));
    else
        d_session = dir(sprintf('/nifvault/procdata/waidmannen/physiology/cells/%s', nameSession));
    end
    setBlockName = {d_session.name}';
    
    if flagSaveFP
        % Save FingerPrinting data
        locFPBlock = find(contains(setBlockName, 'FPrint')>0);
        
        FPData = struct([]);
        for iBlock = 1:length(locFPBlock)
            d_block = dir(fullfile(d_session(locFPBlock(iBlock)).folder, d_session(locFPBlock(iBlock)).name, '*.mat'));
            setFileName = {d_block.name}';
            
            FPData(iBlock).info.name = fullfile(d_session(locFPBlock(iBlock)).folder, d_session(locFPBlock(iBlock)).name);
            
            % load analog data
            load(fullfile(d_block(end).folder, 'tankdat.mat'));
            
            % photodiode & stim ON/OFF
            FPData(iBlock).info.fs_analog = tankdata.streams.Anlg.fs; 
            FPData(iBlock).info.tAxis_ms = [0:1/tankdata.streams.Anlg.fs:(size(tankdata.streams.Anlg.data, 2)-1)/tankdata.streams.Anlg.fs].*1000;
            
            data_pd    = tankdata.streams.Anlg.data(3,:);
            pd_threshold = (max(data_pd)-min(data_pd))/5;
            FPData(iBlock).pd = data_pd;
            FPData(iBlock).pd_threshold = pd_threshold;
            
%             FPData(iBlock).pd_on    = find(diff(data_pd)>pd_threshold)+1;
%             FPData(iBlock).pd_on    = FPData(iBlock).pd_on./FPData(iBlock).info.fs_analog.*1000;
%             FPData(iBlock).pd_off = find(diff(data_pd)<-pd_threshold)+1;
%             FPData(iBlock).pd_off = FPData(iBlock).pd_off./FPData(iBlock).info.fs_analog.*1000;        
            
            % convert the eye signal to visual angles
            clear eye*
            eyedata = tankdata.streams.Anlg.data(1:2,:); % extract eye movement data (in voltage unit)
            eyevector = [eyedata', ones(size(eyedata, 2), 1)]; % transpose and add 1 (constant term) for calculation of matrix product
            tmat = eTform.tdata.T; % affine transformation matrix
            eyedeg = eyevector * tmat; % affine transformation            
            
            FPData(iBlock).eyedeg = eyedeg(:, 1:2);
            
            % Eventcode & timing
            FPData(iBlock).event_code = tankdata.scalars.Evnt.data;
            FPData(iBlock).event_time_ms = tankdata.scalars.Evnt.ts.*1000; % to millisecond
            
            % load cell data
            for iCell = 1:length(setFileName)-1
                
                load(fullfile(d_block(iCell).folder, d_block(iCell).name))
                
                FPData(iBlock).info.cells(iCell) = celldata.info;
                FPData(iBlock).cell_ts_org{iCell, 1} = celldata.ts; %
%                 FPData(iBlock).cell_ts_aligned{iCell, 1} = celldata.ts - t_movieOn_ms; % align it to the start of the movie
%                 movieData(iBlock).info.t_movieOn_ms = t_movieOn_ms;
                
                fprintf(1, 'Session: %s, Block %d/%d: %s, Cell %d/%d processed \n', nameSession, iBlock, length(locFPBlock), ...
                    d_session(locFPBlock(iBlock)).name, iCell, length(setFileName)-1);
            end
            
        end
        fprintf(1, 'Session: %s ....saving... \n', nameSession)
        switch lower(char(regexp(nameSession, '\D*', 'match')))
            case 'spice'
                nameArea = 'AF';
            case 'wasabi'
                nameArea = 'AM';
            case 'davida'
                nameArea = 'ML';
            case 'mochi'
                nameArea = 'AFaAM';
        end
%         nameArea = 'AM';
        saveFileName = sprintf('/nifvault/procdata/parksh/_macaque/_Marta/%s_%s_FPrint.mat', nameArea, nameSession);
        save(saveFileName, 'FPData')
        fprintf(1, 'Session: %s ....DONE! Data saved as %s \n', nameSession, saveFileName)
    end
    
    
    if flagSaveMovie
        % Save movie data
        locMovieBlock = find(contains(setBlockName, 'Movie')>0);
        
        movieData = struct([]);
        for iBlock = 1:length(locMovieBlock)
            d_block = dir(fullfile(d_session(locMovieBlock(iBlock)).folder, d_session(locMovieBlock(iBlock)).name, '*.mat'));
            setFileName = {d_block.name}';
            
            % load analog data
            load(fullfile(d_block(end).folder, 'tankdat.mat'));
            
            % get the timestamp for the start of the movie
            t_movieOn_ms = tankdata.scalars.Evnt.ts(find(tankdata.scalars.Evnt.data == 100, 1)).*1000; % convert timestamp for event code 100 to millisecond
            
            % convert the eye signal to visual angles
            clear eye*
            eyedata = tankdata.streams.Anlg.data(1:2,:); % extract eye movement data (in voltage unit)
            eyevector = [eyedata', ones(size(eyedata, 2), 1)]; % transpose and add 1 (constant term) for calculation of matrix product
            tmat = eTform.tdata.T; % affine transformation matrix
            eyedeg = eyevector * tmat; % affine transformation
            
            
            movieData(iBlock).eyedeg = eyedeg(:, 1:2);
            movieData(iBlock).info.eye.fs = tankdata.streams.Anlg.fs;
            movieData(iBlock).info.eye.tAxis_ms = [0:1/tankdata.streams.Anlg.fs:(length(eyedeg)-1)/tankdata.streams.Anlg.fs].*1000;
            movieData(iBlock).info.eye.tAxis_ms_aligned = movieData(iBlock).info.eye.tAxis_ms - t_movieOn_ms;
            
            % load cell data
            for iCell = 1:length(setFileName)-1
                
                load(fullfile(d_block(iCell).folder, d_block(iCell).name))
                
                movieData(iBlock).info.cells(iCell) = celldata.info;
                movieData(iBlock).cell_ts_org{iCell, 1} = celldata.ts; %
                movieData(iBlock).cell_ts_aligned{iCell, 1} = celldata.ts - t_movieOn_ms; % align it to the start of the movie
                movieData(iBlock).info.t_movieOn_ms = t_movieOn_ms;
                
                fprintf(1, 'Session: %s, Block %d/%d: %s, Cell %d/%d processed \n', nameSession, iBlock, length(locMovieBlock), ...
                    d_session(locMovieBlock(iBlock)).name, iCell, length(setFileName)-1);
            end
            
        end
        fprintf(1, 'Session: %s ....saving... \n', nameSession)
        switch lower(char(regexp(nameSession, '\D*', 'match')))
            case 'spice'
                nameArea = 'AF';
            case 'wasabi'
                nameArea = 'AM';
            case 'davida'
                nameArea = 'ML';
            case 'mochi'
                nameArea = 'AFaAM';
        end
        saveFileName = sprintf('/nifvault/procdata/parksh/_macaque/_Marta/%s_%s.mat', nameArea, nameSession);
        save(saveFileName, 'movieData')
        fprintf(1, 'Session: %s ....DONE! Data saved as %s \n', nameSession, saveFileName)
    end
    
end




% load /procdata/koyanok/physiology/cells/Spice180120/Movie10_4/tankdat.mat; % load analog data
%
% eyedata = tankdata.streams.Anlg.data(1:2,:); % extract eye movement data (in voltage unit)
% eyevector = [eyedata’, ones(size(eyedat, 1))]; % transpose and add 1 (constant term) for calculation of matrix product
% tmat = eTform.tdata.T; % affine transformation matrix
% eyedeg = eyevector * tmat; % affine transformation
% plot(eyedeg(:,1), eyedeg(:,2)); % plot eye movement data