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
setSession = {'Spice180124', 'Spice180126'}; %{'Wasabi190304', 'Wasabi190305', 'Wasabi190306'}; %{'Davida180723', 'Davida180724', 'Davida180725'}; %{'Spice180124', 'Spice180126'};

for iSession = 1:length(setSession)
    nameSession = setSession{iSession};
    typeTankData = 1;
    dateSession = nameSession(end-5:end);
    if str2num(dateSession) > 190700
        typeTankData = 2;
    end
    
    d_session = dir(sprintf('/procdata/koyanok/physiology/cells/%s', nameSession));
%     d_session = dir(sprintf('/procdata/waidmannen/physiology/cells/%s', nameSession));
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
        
        
        movieData(iBlock).eyedeg = eyedeg(:, 1:2);
        movieData(iBlock).pupildata = pupildata;
        movieData(iBlock).info.eye.fs = tankdata.streams.Anlg.fs;
        movieData(iBlock).info.eye.tAxis_ms = [0:1/tankdata.streams.Anlg.fs:(length(eyedeg)-1)/tankdata.streams.Anlg.fs].*1000;
        movieData(iBlock).info.eye.tAxis_ms_aligned = movieData(iBlock).info.eye.tAxis_ms - t_movieOn_ms;
        movieData(iBlock).info.t_movieOn_ms = t_movieOn_ms;
        movieData(iBlock).info.nameBlock = fullfile(d_session(locMovieBlock(iBlock)).folder, d_session(locMovieBlock(iBlock)).name);
        movieData(iBlock).info.movieID = d_session(locMovieBlock(iBlock)).name(end);
        
        % quick figure
        durMovie_ms = 300000;
        figure;
        set(gcf, 'Color', 'w')
        plot(movieData(iBlock).info.eye.tAxis_ms, movieData(iBlock).pupildata)
        hold on
        line([movieData(iBlock).info.t_movieOn_ms movieData(iBlock).info.t_movieOn_ms], get(gca, 'YLim'), 'Color', 'm')
        line([movieData(iBlock).info.t_movieOn_ms movieData(iBlock).info.t_movieOn_ms]+durMovie_ms, get(gca, 'YLim'), 'Color', 'm')
        ylim([-8 -5])        
        ylabel('Pupil size')
        xlabel('Time (ms)')
        title(sprintf('%s: Movie %s', nameSession, d_session(locMovieBlock(iBlock)).name(end)))        
    end
    
end

rs = resample(double(pupildata), 10000000000, 30517578125);
















