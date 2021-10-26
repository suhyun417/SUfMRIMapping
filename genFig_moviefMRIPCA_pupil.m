% genFig_moviefMRIPCA_pupil.m
%
% 2021/07/12 SHP: investigate whether pupil size changes are related to
% movie fMRI PCs

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

%% Load data
load('/procdata/parksh/_macaque/ArtAva_movieTS_fMRI_concatTS_pca.mat', 'S'); % fMRI TS & PCs
% load('/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat', 'infoTS_subj', 'matTS_FP'); % neuronal time courses
%
% % face selectivity
% load('/procdata/parksh/_macaque/multipleFP_fsi.mat')


%% re-form principal components
for iPC = 1:3
    for iS = 1:2
        
        tt = reshape(S(iS).resultsPCA_concat_brainmask.coeff(:,iPC), 118, 3);
        tt = cat(1, NaN(7, 3), tt);
        
        PC{iPC}(:,iS) = tt(:); %S(iS).resultsPCA_concat_brainmask.coeff(:,iPC);
    end
end

figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iPC = 1:3
    SP(iPC) = subplot(3,1,iPC);
    P{iPC} = plot(2.4:2.4:900, PC{iPC}, '-', 'LineWidth', 2);
    axis tight
end

%% Load pupil size change from one animal
% load transformation matrix
load /projects/koyanok/Codes/BehaviorControl/MonkeyLogic/Volt2DegRig1.mat;
setSession{1} = {'Spice180124', 'Spice180126'};
setSession{2} = {'Wasabi190304', 'Wasabi190305', 'Wasabi190306'};
setSession{3} = {'Davida180723', 'Davida180724', 'Davida180725'};
S = struct([]);

for iSubj = 1:length(setSession)
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
        if iSubj > 1
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

% save file
saveFileName = 'eyeSignal_Movie123_ePhys_SpiWasDav.mat';
save(fullfile('/procdata/parksh/_macaque/', saveFileName), 'S', '-v7.3')

%
matPupilData_subj = [];
for idMovie = 1:3
for iSubj = 1:3
    curM = S(iSubj).Session(1).movieData;
    [sortedMovieIDs, indBlock] = sort(str2num(cat(1, curM.movieID)));
    locCurMovieBlock = indBlock(sortedMovieIDs==idMovie);
    matPupilData_subj = cat(2, matPupilData_subj, curM(locCurMovieBlock).pupildata_ms_movie);
end
end
matR_subj = corr(matPupilData_subj, 'rows', 'complete', 'type', 'Spearman');


[sortedMovieIDs, indBlock] = sort(str2num(cat(1, movieData.movieID)));
idMovie = 1;
locCurMovieBlock = indBlock(sortedMovieIDs==idMovie);

matPupilData = cat(2, movieData(locCurMovieBlock).pupildata_ms_movie);
meanPupilData = nanmean(matPupilData, 2);

matR = corr(matPupilData, 'rows', 'complete', 'type', 'Spearman');

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(matR)
set(gca, 'CLim', [-1 1])
% iBlock = 1;
% idBlock = locMovieBlock(iBlock);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(matPupilData)
hold on
plot(meanPupilData, 'k-', 'LineWidth', 2)

%% pupil regressor
% MION function
addpath('/library/matlab_utils');
TR=2.4;
k = gampdf([-40:TR:40],4,2);

pupildata_TR = resample(meanPupilData, 1, 2400);

matNeuralRGR = NaN(nt, length(validC));
for iChan = 1:length(validC)
    neuralrgrs=[];
    for iMov = 1:length(indMovieNeuron)
        curNeuralTC = S(validC(iChan), indMovieNeuron(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
        curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
        curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
        curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
        neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
    end
    matNeuralRGR(:,iChan) = neuralrgrs';
end












