% saveMovieData_SU_EyeSignal_ForHarish.m
% 2022/10/05 SHP
%
%   - This script needs customization mainly for retrieving directory names
%   - This structure is only true for the cells with face selectivity
%   measures (image responses) except for the cells from Mat (monkey
%   Matcha): this is Kenji's data structure and the other cells were recorded by other people
%   - To gather eye data corresponding to a particular trial of a
%   particular cell, start from
%       : load('/nifvault/procdata/parksh/_macaque/_Harish/TS_movie123_30fps.mat')
%   - Then using parameters saved in matSDF to load the corresponding file
%   containing the information for the eye data directory in this way
%       : load(sprintf('/nifvault/procdata/parksh/_macaque/%s/%s', ...
%       matSDF(iSubj).nameSubj, matSDF(iSubj).FR_30fps(iCell,iMovie).SUfilename))
%   - Use dat.tank and dat.trial to find the corresponding eye data
%   - below is what is in the analog data
% tankdata.streams.Anlg.data : x , y, photodiode, reward, pupil diameter
%                                                   (from channel 1 to 5)
% tankdata.streams.Anlg.fs : analog signal sampling rate
% tankdata.scalars.Evnt.data : event code (100 is the start of the movie)
% tankdata.scalars.Evnt.ts: timestamp of event code. in seconds.
% celldata.ts: timestamp of spikes. in milliseconds.
% celldata.ts and tankdata.streams.Anlg, tankdata.scalars.Evnt.ts are
% aligned to the start of the recording



clear all;


% load transformation matrix
load /projects/koyanok/Codes/BehaviorControl/MonkeyLogic/Volt2DegRig1.mat;
%%% NOTE: throughout the code you need to add "nifvault" at the beginning

setSession{1} = {'Spice180124', 'Spice180126'};
setSession{2} = {'Wasabi190304', 'Wasabi190305', 'Wasabi190306'};
setSession{3} = {'Davida180723', 'Davida180724', 'Davida180725'};
S = struct([]);

for iSubj = 1:length(setSession) %%% NOTE: this "iSubj" is different from the matSDF subject indices
    curSetSession = setSession{iSubj}; %{'Spice180124', 'Spice180126'}; %{'Wasabi190304', 'Wasabi190305', 'Wasabi190306'}; %{'Davida180723', 'Davida180724', 'Davida180725'}; %{'Spice180124', 'Spice180126'};
    
    Session = struct([]);
    for iSession = 1:length(curSetSession)
        nameSession = curSetSession{iSession};
        typeTankData = 1;
        dateSession = nameSession(end-5:end);
        if str2num(dateSession) > 190700
            typeTankData = 2;
        end
        
        fprintf(1, 'Processing Subject #%d/%d: Session #%d/%d (%s) ....... \n', iSubj, length(setSession), iSession, length(curSetSession), ...
            nameSession);
        
        d_session = dir(sprintf('/procdata/koyanok/physiology/cells/%s', nameSession));
        if iSubj > 1 %%% NOTE: this is because Wasabi and Davida's files are processed by Elena Waidmann
        d_session = dir(sprintf('/procdata/waidmannen/physiology/cells/%s', nameSession));
        end
        setBlockName = {d_session.name}';
        locMovieBlock = find(contains(setBlockName, 'Movie')>0);
        %%% NOTE: above two lines to get ALL the Movie trials available from
        %%% this daily session. You can use "dat.trial" field to get specific block 
        %%% that corresponds to your cell's particular trial
        
        movieData = struct([]);
        for iBlock = 1:length(locMovieBlock)
            d_block = dir(fullfile(d_session(locMovieBlock(iBlock)).folder, d_session(locMovieBlock(iBlock)).name, '*.mat'));
            setFileName = {d_block.name}';
            
            % load analog data
            load(fullfile(d_block(end).folder, 'tankdat.mat'));
            
            switch typeTankData
                case 1
                    pupildata = tankdata.streams.Anlg.data(5, :)';
                case 2
                    pupildata = tankdata.streams.Anlg2.data(1, :)';
            end
            
            % get the timestamp for the start of the movie
            t_movieOn_ms = tankdata.scalars.Evnt.ts(find(tankdata.scalars.Evnt.data == 100, 1)).*1000; % convert timestamp for event code 100 to millisecond
            
            % convert the eye signal to visual angles
            clear eye*
            eyedata = tankdata.streams.Anlg.data(1:2,:); % extract eye movement data (in voltage unit)
            eyevector = [eyedata', ones(size(eyedata, 2), 1)]; % transpose and add 1 (constant term) for calculation of matrix product
            tmat = eTform.tdata.T; % affine transformation matrix
            eyedeg = eyevector * tmat; % affine transformation
            
            eyedeg_ms = resample(double(eyedeg(:, 1:2)), 10000000000, 30517578125);
            eyedeg_ms_movie = eyedeg_ms(ceil(t_movieOn_ms):ceil(t_movieOn_ms)+300000-1, :);
            pupildata_ms = resample(double(pupildata), 10000000000, 30517578125);
            pupildata_ms_movie = pupildata_ms(ceil(t_movieOn_ms):ceil(t_movieOn_ms)+300000-1);
            %%% NOTE: resample in your desired resolution by changing the
            %%% second input of the "resample" functions above. Currently
            %%% in ms resolution.
            
            % Blink correction
            blinks_data_positions = [];
            sampling_rate_hz = 1000;
            blink_threshold_voltage = -8; % voltage threshold for loss of pupil
            blinks_data_positions = based_noise_blinks_detection_shp(pupildata_ms_movie, sampling_rate_hz, blink_threshold_voltage);	% get blink positions using the noise-based approach
            if(~isempty(blinks_data_positions) && blinks_data_positions(1)==0)
                blinks_data_positions(1) = 1; %[];
            end
            
            blinks_data_positions = reshape(blinks_data_positions, 2, length(blinks_data_positions)./2)';
            for iBlink = 1:size(blinks_data_positions, 1)
                pupildata_ms_movie(blinks_data_positions(iBlink,1)+1:blinks_data_positions(iBlink,2)) = NaN;
                eyedeg_ms_movie(blinks_data_positions(iBlink,1)+1:blinks_data_positions(iBlink,2), :) = NaN;
            end
            %%% NOTE: above is my attempt to detect & remove the blinks. I
            %%% put the function "based_noise_blinks_detection_shp" in
            %%% /_Harish directory too
            
            movieData(iBlock).eyedeg_ms = eyedeg_ms;
            movieData(iBlock).pupildata_ms = pupildata_ms;
            movieData(iBlock).blinks_data_positions = blinks_data_positions;
            movieData(iBlock).eyedeg_ms_movie = eyedeg_ms_movie;
            movieData(iBlock).pupildata_ms_movie = pupildata_ms_movie;
            movieData(iBlock).info.eye.fs = tankdata.streams.Anlg.fs;
            movieData(iBlock).info.eye.tAxis_ms = [0:1/tankdata.streams.Anlg.fs:(length(eyedeg)-1)/tankdata.streams.Anlg.fs].*1000;
            movieData(iBlock).info.eye.tAxis_ms_aligned = movieData(iBlock).info.eye.tAxis_ms - t_movieOn_ms;
            movieData(iBlock).info.t_movieOn_ms = t_movieOn_ms;
            movieData(iBlock).info.nameBlock = fullfile(d_session(locMovieBlock(iBlock)).folder, d_session(locMovieBlock(iBlock)).name);
            movieData(iBlock).movieID = d_session(locMovieBlock(iBlock)).name(end);
            
            %         % quick figure
            %         durMovie_ms = 300000;
            %         figure;
            %         set(gcf, 'Color', 'w')
            %         plot(movieData(iBlock).info.eye.tAxis_ms, movieData(iBlock).pupildata)
            %         hold on
            %         line([movieData(iBlock).info.t_movieOn_ms movieData(iBlock).info.t_movieOn_ms], get(gca, 'YLim'), 'Color', 'm')
            %         line([movieData(iBlock).info.t_movieOn_ms movieData(iBlock).info.t_movieOn_ms]+durMovie_ms, get(gca, 'YLim'), 'Color', 'm')
            %         ylim([-8 -5])
            %         ylabel('Pupil size')
            %         xlabel('Time (ms)')
            %         title(sprintf('%s: Movie %s', nameSession, d_session(locMovieBlock(iBlock)).name(end)))
        end
        
        Session(iSession).name = nameSession;
        Session(iSession).movieData = movieData;
        
    end
    S(iSubj).Session = Session;
    
end



