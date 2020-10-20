% computeCellParam_waveform.m
%
% 2020/10/14 SHP
% This code computes and saves the waveform parameters (waveforms, fitted lines, peak-to-trough widths) 
% for all 401 neurons from AF, AM, AAM, peri-AM.
% The response time course of these neurons to movies are saved in 
% '/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat'

clear all;

%% Set the directory and load a file
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

addpath('/projects/parksh/_toolbox/Waveform_Consistency-V1.1/')

load('/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat', 'infoTS_subj', 'matTS_all')


%%
setMovIDs = [1 2 3];
resultsWF_subj = struct([]);
for iSubj = 1:length(infoTS_subj)
    
    nameSubject = infoTS_subj(iSubj).nameSubj;
    dirDataNeural = fullfile(dirProcdata, '/parksh/_macaque', nameSubject);
    
    listSUchannelID = infoTS_subj(iSubj).validChanID;
    unitchar = struct([]);
    
    for iM = 1:length(setMovIDs)
        movID = setMovIDs(iM); %1;
        d_n = dir(fullfile(dirDataNeural, sprintf('*mov%dsig*.mat', movID)));
        
        listMatSUFile = {d_n.name}';
        
        for iChan = 1:length(listSUchannelID)
            curChanID = ['sig' listSUchannelID{iChan} '.'];
            indFile = find(contains(listMatSUFile, curChanID)>0);
            
            load(fullfile(d_n(indFile).folder, d_n(indFile).name))
            
            if iSubj < 4 % old dataset with different structure
                resultsWF_subj(iSubj).nameSubj = infoTS_subj(iSubj).nameSubj;
                resultsWF_subj(iSubj).validChanID = infoTS_subj(iSubj).validChanID;
                resultsWF_subj(iSubj).meanWaveP2T = NaN(size(resultsWF_subj(iSubj).validChanID, 1), 1);
                resultsWF_subj(iSubj).steWaveP2T = NaN(size(resultsWF_subj(iSubj).validChanID, 1), 1);
                
                %% need to compute number of spikes.. probably trial by trial?? darn
                % Toroid: data.s{iTrial}
                
            elseif iSubj == 5 % Matchat's case: different data structure
                % number of spikes : size(dat.s{iTrial}, 1)
                % no waveform information
            else
            
            
            for iTrial = 1:length(dat.waveform)
                
                if size(dat.waveform{iTrial}, 2)<5
                    unitchar(iChan, iM).numSpikes(iTrial, 1) = size(dat.waveform{iTrial}, 2);
                    unitchar(iChan, iM).matWaveP2T(iTrial, 1) = NaN;
                    unitchar(iChan, iM).fitWaveResult(iTrial,1).xaxis_microsec = NaN;
                    unitchar(iChan, iM).fitWaveResult(iTrial,1).splineWaveform = NaN;
                    continue;
                end
                
                wf = nanmean(dat.waveform{iTrial}');
                [p2t, xaxis, yy] = getWaveformPeakToTrough(wf);
                %             figure(30);
                %             plot(xaxis, yy); hold on
                
                unitchar(iChan, iM).movID = movID;
                unitchar(iChan, iM).numSpikes(iTrial, 1) = size(dat.waveform{iTrial}, 2);
                unitchar(iChan, iM).matWaveP2T(iTrial, 1) = p2t;
                unitchar(iChan, iM).fitWaveResult(iTrial,1).xaxis_microsec = xaxis;
                unitchar(iChan, iM).fitWaveResult(iTrial,1).splineWaveform = yy;
            end
            %         figure(30); clf;
        end
    end
    
    meanWaveP2T = []; steWaveP2T = [];
    meanNumSpikes = []; steNumSpikes = [];
    unitchar_avgmovie = struct([]);
    for iChan = 1:size(unitchar, 1)
        unitchar_avgmovie(iChan, 1).matWaveP2T = cat(1, unitchar(iChan, :).matWaveP2T);
        unitchar_avgmovie(iChan, 1).matNumSpikes = cat(1, unitchar(iChan, :).numSpikes);
        meanWaveP2T(iChan, 1) = nanmean(unitchar_avgmovie(iChan, 1).matWaveP2T);
        steWaveP2T(iChan, 1) = nanstd(unitchar_avgmovie(iChan, 1).matWaveP2T)./sqrt(size(unitchar_avgmovie(iChan, 1).matWaveP2T, 1));
        meanNumSpikes(iChan, 1) = nanmean(unitchar_avgmovie(iChan, 1).matNumSpikes);
        steNumSpikes(iChan, 1) = nanstd(unitchar_avgmovie(iChan, 1).matNumSpikes)./sqrt(size(unitchar_avgmovie(iChan, 1).matNumSpikes, 1));
    end
    
    resultsWF_subj(iSubj).nameSubj = infoTS_subj(iSubj).nameSubj;
    resultsWF_subj(iSubj).validChanID = infoTS_subj(iSubj).validChanID;
    resultsWF_subj(iSubj).unitchar_eachMovie = unitchar;
    resultsWF_subj(iSubj).unitchar_avgmovie = unitchar_avgmovie;
%     resultsWF_subj(iSubj).matWaveP2T = matWaveP2T;
    resultsWF_subj(iSubj).meanWaveP2T = meanWaveP2T; %nanmean(matWaveP2T);
    resultsWF_subj(iSubj).steWaveP2T = steWaveP2T; %nanstd(matWaveP2T)./sqrt(size(matWaveP2T, 1));
%     resultsWF_subj(iSubj).matNumSpikes = matNumSpikes;
    resultsWF_subj(iSubj).meanNumSpikes = meanNumSpikes; %mean(matNumSpikes);
    resultsWF_subj(iSubj).steNumSpikes = steNumSpikes; %std(matWaveP2T)./sqrt(size(matWaveP2T, 1));
    
end

figure;
errorbar(resultsWF_subj(iSubj).meanWaveP2T, resultsWF_subj(iSubj).steWaveP2T)
text(1:length(resultsWF_subj(iSubj).meanWaveP2T), resultsWF_subj(iSubj).meanWaveP2T, listSUchannelID)
set(gca, 'XTick', [])
title(sprintf('%s: spikes during 3 movies', infoTS_subj(iSubj).nameSubj))
ylabel('peak-to-trough width')

figure;
for iM = 1:3
subplot(3,1,iM)
tempMat = cat(2, unitchar(:, iM).matWaveP2T);
plot(nanmean(tempMat), 'o-');
text(1:24, nanmean(tempMat), listSUchannelID);
title(sprintf('%s: movie %d', infoTS_subj(iSubj).nameSubj, iM))
end

%% Load the waveform and characterize one by one 
for iChan = 1:length(listSUchannelID)
  
  cur_cellID = listSUchannelID{iChan};
  indCurChan = find(strcmp(cur_cellID, listSUchannelID)>0);
  
  for iMov = 1:length(setMovIDs)
      cur_movID = setMovIDs(iMov);
      if ~indDataMov(indCurChan, iMov) 
          fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
          fprintf(1, 'Skip to the next movie/cell \n')
          continue; 
      end
      
      % Find a relevant file
      filename = char(listMatSUFile((strcmp(cur_cellID, listSU_all)+strcmp(num2str(cur_movID), listMov_all))==2));
      fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
      fprintf(1, 'filename: %s\n', filename)
      
      load(fullfile(dirDataNeural, filename))
  
      for iTrial = 1:length(dat.waveform)
          
          wf = nanmean(dat.waveform{iTrial}');
          [p2t, xaxis, yy] = getWaveformPeakToTrough(wf);
          %         figure(30);
          %         plot(xaxis, yy); hold on
          
          unitchar(iChan, iMov).matWaveP2T(iTrial, 1) = p2t;
          unitchar(iChan, iMov).fitWaveResult(iTrial,1).xaxis_microsec = xaxis;
          unitchar(iChan, iMov).fitWaveResult(iTrial,1).splineWaveform = yy;
      end
  end
end


%   count = 0;
  for iMov = 1:length(setMovIDs) % should concatenate all those movies in order
      
      cur_movID = setMovIDs(iMov);
      
      if ~indDataMov(indCurChan, iMov), 
          fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
          fprintf(1, 'Skip to the next movie/cell \n')
          continue; 
      end
      
      % Find a relevant file
      filename = char(listMatSUFile((strcmp(cur_cellID, listSU_all)+strcmp(num2str(cur_movID), listMov_all))==2));
      fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
      fprintf(1, 'filename: %s\n', filename)

      FR_dT(iChan, iMov).cellID = cur_cellID;
      FR_dT(iChan, iMov).movID = cur_movID;
      FR_dT(iChan, iMov).SUfilename = filename;
      [FR_dT(iChan, iMov).matFR{1}, FR_dT(iChan, iMov).mnFR] = computeMeanFR(dirDataNeural, filename, sizeTimeBin_sec); % averaged SDF across trials

%       [S(iChan, iMov).matsdf{1}, S(iChan, iMov).mnsdf] = computeMeanSDF(dirDataNeural, filename); % averaged SDF across trials
%       [S(iChan, iMov).rgrsMION, S(iChan, iMov).rgrsMION_resample] = computeMIONrgrs(SDF_dT(iChan, iMov).mnsdf, timeResNeural_sec, fMRI_TR_sec);

      
  end
      
  
% end

% [p2t, xaxis, yy] = getWaveformPeakToTrough(wf) 

% %% Idnetify available units
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
