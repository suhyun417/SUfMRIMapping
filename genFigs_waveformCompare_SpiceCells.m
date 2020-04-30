

%% Set the directory
ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS';
    dirProcdata = '/Volumes/PROCDATA';
    dirLibrary = '/Volumes/LIBRARY';
else % on virtual machine
    dirProjects = '/projects';
    dirProcdata = '/procdata';
    dirLibrary = '/library';
end

nameSubject = 'Spi';
dirDataNeural = fullfile(dirProcdata, '/parksh', nameSubject);

setMovIDs = [1 2 3];



% ID of cells from Spice FingerPrinting data on 11/28/16
faceCells = {'005_1' '009_1' '014_1' '015_1' '017_1' '020_2' '024_1' '025_1' '034_1' '036_1'};
facebodyCells = {'001_1' '002_1' '003_1' '008_1' '037_2' '046_2' '049_1' '060_2' '061_1' '062_1'};


% first, start from the same channel from both dataset
% start from the finger printing dataset, one-by-one
% see whether movie data has a unit from the same channel
% then pulled out two dataset
% first waveform
% then ISI
length(faceCells)

for iC=1:length(faceCells)
    idCurCell=char(faceCells(iC));
    % curChan = str2num(idCurCell(1:3));
    
    cellCompare_face(iC).cellID = faceCells{iC};
    
    % finger printing data from 11/28/16
    load(sprintf('/procdata/koyanok/physiology/cells/Spice161128/FPrint1/%s.mat', idCurCell))
    isExist(1,1) = 1;
    ISI{1} = diff(celldata.ts);
    wf{1} = nanmean(celldata.wf,2);
    
    % movie data from 11/20/16: check possible two units
    if ~exist(sprintf('/procdata/koyanok/physiology/cells/Spice161120/1_mov1_1/%s_1.mat', idCurCell(1:3)))
        isExist(2,1) = 0;
        ISI{2} = [];
        wf{2} = [];
    else
        load(sprintf('/procdata/koyanok/physiology/cells/Spice161120/1_mov1_1/%s_1.mat', idCurCell(1:3)))
        isExist(2,1) = 1;
        ISI{2} = diff(celldata.ts);
        wf{2} = nanmean(celldata.wf,2);
    end
    if ~exist(sprintf('/procdata/koyanok/physiology/cells/Spice161120/1_mov1_1/%s_2.mat', idCurCell(1:3)))
        isExist(3,1) = 0;
        ISI{3} = [];
        wf{3} = [];
    else
        load(sprintf('/procdata/koyanok/physiology/cells/Spice161120/1_mov1_1/%s_2.mat', idCurCell(1:3)))
        isExist(3,1) = 1;
        ISI{3} = diff(celldata.ts);
        wf{3} = nanmean(celldata.wf,2);
    end
    
    if sum(isExist)>1
        matR = corrcoef(cell2mat(wf));
    else
        matR = [];
    end
    
    cellCompare_face(iC).isExist = isExist;
    cellCompare_face(iC).ISI = ISI;
    cellCompare_face(iC).wf = wf;
    cellCompare_face(iC).matR = matR;

end


for iC=1:length(facebodyCells)
    idCurCell=char(facebodyCells(iC));
    % curChan = str2num(idCurCell(1:3));
    
    cellCompare_facebody(iC).cellID = facebodyCells{iC};
    
    % finger printing data from 11/28/16
    load(sprintf('/procdata/koyanok/physiology/cells/Spice161128/FPrint1/%s.mat', idCurCell))
    isExist(1,1) = 1;
    ISI{1} = diff(celldata.ts);
    wf{1} = nanmean(celldata.wf,2);
    
    % movie data from 11/20/16: check possible two units
    if ~exist(sprintf('/procdata/koyanok/physiology/cells/Spice161120/1_mov1_1/%s_1.mat', idCurCell(1:3)))
        isExist(2,1) = 0;
        ISI{2} = [];
        wf{2} = [];
    else
        load(sprintf('/procdata/koyanok/physiology/cells/Spice161120/1_mov1_1/%s_1.mat', idCurCell(1:3)))
        isExist(2,1) = 1;
        ISI{2} = diff(celldata.ts);
        wf{2} = nanmean(celldata.wf,2);
    end
    if ~exist(sprintf('/procdata/koyanok/physiology/cells/Spice161120/1_mov1_1/%s_2.mat', idCurCell(1:3)))
        isExist(3,1) = 0;
        ISI{3} = [];
        wf{3} = [];
    else
        load(sprintf('/procdata/koyanok/physiology/cells/Spice161120/1_mov1_1/%s_2.mat', idCurCell(1:3)))
        isExist(3,1) = 1;
        ISI{3} = diff(celldata.ts);
        wf{3} = nanmean(celldata.wf,2);
    end
    
    if sum(isExist)>1
        matR = corrcoef(cell2mat(wf));
    else
        matR = [];
    end
    
    cellCompare_facebody(iC).isExist = isExist;
    cellCompare_facebody(iC).ISI = ISI;
    cellCompare_facebody(iC).wf = wf;
    cellCompare_facebody(iC).matR = matR;

end

catSetR=[]; setCellID=[];
for iC = 1:length(faceCells)
    if ~isempty(cellCompare_face(iC).matR)
        setR = cellCompare_face(iC).matR(2:2*length(cellCompare_face(iC).matR)-2);
        catSetR = cat(1, catSetR, max(setR));
        setCellID = cat(1, setCellID, cellCompare_face(iC).cellID);
    else
        catSetR = cat(1, catSetR, 0);
        setCellID = cat(1, setCellID, cellCompare_face(iC).cellID);
    end
end

waveformCompare(1).info = 'face cells';
waveformCompare(1).setMaxR = catSetR;
waveformCompare(1).setCellID = setCellID;

catSetR=[]; setCellID=[];
for iC = 1:length(facebodyCells)
    if ~isempty(cellCompare_facebody(iC).matR)
        setR = cellCompare_facebody(iC).matR(2:2*length(cellCompare_facebody(iC).matR)-2);
        catSetR = cat(1, catSetR, max(setR));
        setCellID = cat(1, setCellID, cellCompare_facebody(iC).cellID);
    else
        catSetR = cat(1, catSetR, 0);
        setCellID = cat(1, setCellID, cellCompare_face(iC).cellID);
    end
end

waveformCompare(2).info = 'face & body cells';
waveformCompare(2).setMaxR = catSetR;
waveformCompare(2).setCellID = setCellID;

dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');
threshR = 0.99;
indFaceCells = find(waveformCompare(1).setMaxR>threshR);
for iFaceCell = 1:length(indFaceCells)
    curInd = indFaceCells(iFaceCell);
    
    figure;
    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [200 200 1000 410])
    subplot(1, 2, 1);
    plot(cell2mat(cellCompare_face(curInd).wf), 'LineWidth', 2)
    box off
    subplot(1,2,2);
    plot(zscore(cell2mat(cellCompare_face(curInd).wf)), 'LineWidth', 2)
    box off
    print(gcf, fullfile(dirFig, sprintf('Revision2_SpiceCells_waveform_%s', cellCompare_face(curInd).cellID)), '-depsc');
end

indFaceBodyCells = find(waveformCompare(2).setMaxR>threshR);
for iFaceBodyCell = 1:length(indFaceBodyCells)
    curInd = indFaceBodyCells(iFaceBodyCell);
    
    figure;
    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [200 200 1000 410])
    subplot(1, 2, 1);
    plot(cell2mat(cellCompare_facebody(curInd).wf), 'LineWidth', 2)
    box off
    subplot(1,2,2);
    plot(zscore(cell2mat(cellCompare_facebody(curInd).wf)), 'LineWidth', 2)
    box off
    print(gcf, fullfile(dirFig, sprintf('Revision2_SpiceCells_waveform_%s', cellCompare_facebody(curInd).cellID)), '-depsc');
end


% %% Idnetify available units from movie data on 11/20/16 and 11/21/16
% d_n = dir(fullfile(dirDataNeural, '*sig*.mat'));
% [listMatSUFile{1:length(d_n), 1}] = deal(d_n.name);
% 
% % list of channel IDs and movie IDs of each file
% listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')'; % list of channels
% listMov_all = regexp(cat(2,d_n.name), '\d*(?=sig)', 'match')'; % list of movies
% 
% % Since not every channel has data for every movie, 
% % the "data presence matrix (Channels x MovieIDs)" is generated here
% % for which channel has which movie's data: 1 for data, 0 for no data
% listSUchannelID = unique(listSU_all);
% 
% indDataMov=[];
% for iCh=1:length(listSUchannelID) % go throucgh channel-by-channel     
%      tempListMov{iCh} = regexp([listMatSUFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
%      % Get indices for common movies across cells
%      indDataMov(iCh, :) = ismember(setMovIDs, str2num(char(tempListMov{iCh}))'); % data presence matrix (Channels x MovieIDs) 
% end
% 
% %% Load the waveform and characterize one by one 
% for iChan = 1:length(listSUchannelID)
%   
%   cur_cellID = listSUchannelID{iChan};
%   indCurChan = find(strcmp(cur_cellID, listSUchannelID)>0);
%   
%   for iMov = 1:length(setMovIDs)
%       cur_movID = setMovIDs(iMov);
%       if ~indDataMov(indCurChan, iMov), 
%           fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
%           fprintf(1, 'Skip to the next movie/cell \n')
%           continue; 
%       end
%       
%       % Find a relevant file
%       filename = char(listMatSUFile((strcmp(cur_cellID, listSU_all)+strcmp(num2str(cur_movID), listMov_all))==2));
%       fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
%       fprintf(1, 'filename: %s\n', filename)
%       
%       load(fullfile(dirDataNeural, filename))
%   
%       for iTrial = 1:length(dat.waveform)
%           
%           wf = nanmean(dat.waveform{iTrial}');
%           [p2t, xaxis, yy] = getWaveformPeakToTrough(wf);
%           %         figure(30);
%           %         plot(xaxis, yy); hold on
%           
%           unitchar(iChan, iMov).matWaveP2T(iTrial, 1) = p2t;
%           unitchar(iChan, iMov).fitWaveResult(iTrial,1).xaxis_microsec = xaxis;
%           unitchar(iChan, iMov).fitWaveResult(iTrial,1).splineWaveform = yy;
%       end
%   end
% end
% 
% 
% %   count = 0;
%   for iMov = 1:length(setMovIDs) % should concatenate all those movies in order
%       
%       cur_movID = setMovIDs(iMov);
%       
%       if ~indDataMov(indCurChan, iMov), 
%           fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
%           fprintf(1, 'Skip to the next movie/cell \n')
%           continue; 
%       end
%       
%       % Find a relevant file
%       filename = char(listMatSUFile((strcmp(cur_cellID, listSU_all)+strcmp(num2str(cur_movID), listMov_all))==2));
%       fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
%       fprintf(1, 'filename: %s\n', filename)
% 
%       FR_dT(iChan, iMov).cellID = cur_cellID;
%       FR_dT(iChan, iMov).movID = cur_movID;
%       FR_dT(iChan, iMov).SUfilename = filename;
%       [FR_dT(iChan, iMov).matFR{1}, FR_dT(iChan, iMov).mnFR] = computeMeanFR(dirDataNeural, filename, sizeTimeBin_sec); % averaged SDF across trials
% 
% %       [S(iChan, iMov).matsdf{1}, S(iChan, iMov).mnsdf] = computeMeanSDF(dirDataNeural, filename); % averaged SDF across trials
% %       [S(iChan, iMov).rgrsMION, S(iChan, iMov).rgrsMION_resample] = computeMIONrgrs(SDF_dT(iChan, iMov).mnsdf, timeResNeural_sec, fMRI_TR_sec);
% 
%       
%   end
%       
%   
% end
% 
% [p2t, xaxis, yy] = getWaveformPeakToTrough(wf) 