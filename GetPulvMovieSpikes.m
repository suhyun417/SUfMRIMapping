
%============== Get movie spike data from face selective pulvinar cells

Subject     = 'Dexter';
DataDir     = '/procdata/murphya/Physio/PulvinarFace';
SpikeDir    = ['/procdata/murphya/Physio/MovieData/', Subject];

load(fullfile(DataDir, 'CorrectedCoords_Dexter.mat'));
load(fullfile(DataDir, 'StatsData.mat'));

CellIndx = Stats(1).FaceCellsList; % 1st is Dexter, 2nd is Layla, then Spice, Ava, Matcha
countValidCell = 0;
for c = 1:numel(CellIndx)
    if length(Data) < 1000
        load(fullfile(DataDir, 'CorrectedCoords_Dexter.mat'));
    end

    fprintf(1, '  Processing Cell Index %d, %d/%d\n', CellIndx(c), c, numel(CellIndx))
    Session     = Data(CellIndx(c)).Session;
    Channel     = Data(CellIndx(c)).Channel;
    CellNo      = Data(CellIndx(c)).Cell;
    
    % Load the spike data for this session
    SessionFile = fullfile(SpikeDir, sprintf('MovRep_%s.mat', Session));
    if ~exist(SessionFile, 'file')
        fprintf('Warning! No file was found for %s, channel %d, cell %d!\n', Session, Channel, CellNo);
        continue;
    end
    load(SessionFile);
    
    % Find data from the right channel
    ChIndx = find([Cells.Channel]==Channel & [Cells.CellNo]==CellNo);
    if isempty(ChIndx)
        fprintf('Warning! No movie data was found for %s, channel %d, cell %d!\n', Session, Channel, CellNo);
    else
        countValidCell = countValidCell +1;
        SpikeData(countValidCell).SpikeTimes = Cells(ChIndx).SpikeTimes; %Cells(ChIndx);
        SpikeData(countValidCell).info.Subject = Cells(ChIndx).Subject; 
        SpikeData(countValidCell).info.Session = Cells(ChIndx).Session; 
        SpikeData(countValidCell).info.Channel = Cells(ChIndx).Channel; 
        SpikeData(countValidCell).info.CellNo = Cells(ChIndx).CellNo; 
        SpikeData(countValidCell).info.MovieIDs = Cells(ChIndx).MovieIDs; 
    end
end

save('/procdata/parksh/Dex/_orgData/FaceCellSpike.mat', 'SpikeData')
fprintf(1, '\n Data is saved. \n')
