% tankdata.streams.Anlg.data : x , y, photodiode, reward, pupil diameter
%                                                   (from channel 1 to 5)
% tankdata.streams.Anlg.fs : analog signal sampling rate
% tankdata.scalars.Evnt.data : event code (100 is the start of the movie)
% tankdata.scalars.Evnt.ts: timestamp of event code. in seconds.
% celldata.ts: timestamp of spikes. in milliseconds.
% celldata.ts and tankdata.streams.Anlg, tankdata.scalars.Evnt.ts are
% aligned to the start of the recording

% cd /procdata/parksh/_macaque/_Reza

clear all;

% load transformation matrix
load /projects/koyanok/Codes/BehaviorControl/MonkeyLogic/Volt2DegRig1.mat;

setSession = {'Wasabi190304', 'Wasabi190305', 'Wasabi190306'}; %{'Davida180723', 'Davida180724', 'Davida180725'}; %{'Spice180124', 'Spice180126'};

for iSession = 1:length(setSession)
    nameSession = setSession{iSession};
    
    if contains(nameSession, 'Spice')
        d_session = dir(sprintf('/procdata/koyanok/physiology/cells/%s', nameSession));
    else
        d_session = dir(sprintf('/procdata/waidmannen/physiology/cells/%s', nameSession));
    end
    setBlockName = {d_session.name}';
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
    nameArea = 'AM'; 
    saveFileName = sprintf('/procdata/parksh/_macaque/_Reza/%s_%s.mat', nameArea, nameSession);
    save(saveFileName, 'movieData')
    fprintf(1, 'Session: %s ....DONE! Data saved as %s \n', nameSession, saveFileName)
    
end




% load /procdata/koyanok/physiology/cells/Spice180120/Movie10_4/tankdat.mat; % load analog data
%
% eyedata = tankdata.streams.Anlg.data(1:2,:); % extract eye movement data (in voltage unit)
% eyevector = [eyedata’, ones(size(eyedat, 1))]; % transpose and add 1 (constant term) for calculation of matrix product
% tmat = eTform.tdata.T; % affine transformation matrix
% eyedeg = eyevector * tmat; % affine transformation
% plot(eyedeg(:,1), eyedeg(:,2)); % plot eye movement data