% getRawfilename_movies

flagLocal = 1; %1; % 1 for desktop
switch flagLocal
    case 1 % desktop
        dirData = '/Volumes/PROCDATA/parksh/';
        dirProject = '/Volumes/PROJECTS/parksh/NeuralBOLD/analysis/';
        dirFig = '/Volumes/PROJECTS/parksh/NeuralBOLD/_labNote/_figs/';
    case 0
        dirData = '/procdata/parksh/';
        dirProject = '/projects/parksh/NeuralBOLD/analysis/';
        dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
end

nameSubject = 'Tor';
dirDataNeural = fullfile(dirData, nameSubject);

[DataMatrixMovie, listSUchannelID, listMovie] = genDataMatrix_SU(nameSubject, 0); % Channel x Movie % [indDataMovie, listSUchannelID, listMovie] = genDataMatrix_SU(nameSubj, flagFig)

% % get the list of channels
% d = dir([dirDataNeural, '/*sig*.mat']);
% % tempFileName={};
% % [tempFileName{1:length(d),1}] = deal(d.name);
% listSU_all = regexp(cat(2,d.name), '(?<=sig)\w*', 'match')'; % list of channels
% listMov_all = regexp(cat(2,d.name), '\d*(?=sig)', 'match')'; % list of movies
% 
% listSUchannelID = unique(listSU_all);
listRawFile = {'20120617b', '20120618b', '20120619a', '20120621b', '20120622a', '20120623a', '20120624a'}';

for iChan = 1:length(listSUchannelID)
    idChan = listSUchannelID{iChan};
    matRawIndex = nan(length(listRawFile), length(listMovie));
%     setMovie = listMovie(find(DataMatrixMovie(iChan,:)>0));
    for iMov = 1:length(listMovie)
        idMov = setMovie(iMov);
        if DataMatrixMovie(iChan, iMov)
            movieMatFileName = [dirDataNeural, '/', lower(nameSubject(1)), sprintf('mov%dsig%s.mat', idMov, idChan)];            
            load(movieMatFileName)
            tempRawFilenames = regexp(dat.h.catFnames, '(?<=toroid)\w*', 'match')';
            matRawIndex(:,iMov) = ismember(listRawFile, tempRawFilenames)';
        else
            continue;
        end
    end
    setMatRawIndex{iChan} = matRawIndex;
    
end

figure;
imagesc(setMatRawIndex{1})
colormap(gray(256))
set(gca, 'YTickLabel', listRawFile)
set(gca, 'FontSize', 15)
set(gca, 'XTickLabel', num2str(listMovie))
xlabel('Movie #')
ylabel('Raw file name')
set(gca, 'XAxisLocation', 'top')
set(gcf, 'Color', 'w')
set(gcf, 'PaperPositionMode', 'auto')
% print(gcf, fullfile(dirFig, 'tableMovieRawdatafiles_toroid065a'), '-depsc')
        

% catRawFiles={};
% for iFile = 1:length(movieMatFileName)
% load(fullfile(dirDataNeural, movieMatFileName{iFile}))
% tempRawFilenames = regexp(dat.h.catFnames, '(?<=toroid)\w*', 'match')';
% catRawFiles = cat(1, catRawFiles, tempRawFilenames);
% end
% unique(catRawFiles)