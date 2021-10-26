nameSubjNeural = 'Tor'; % 'Moc'; %'Spi'; %'Was'; % 'Dan'; %'Dav'; %'Was'; %'Ava'; %'Mat'; %'Sig'; %'Rho'; %'Spi';
setMovie = [1 2 3]; % [4 5 6]; %
dirDataHome = '/procdata/parksh/_macaque';
[indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
validC = find(indDataMat*ismember(movieID, setMovie)>0);
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);

load(fullfile(dirDataNeural, [nameSubjNeural, '_movieTS_SU_indMov.mat']), 'paramSDF')
validCellID=paramSDF.setCellIDs(validC);
cellNum=mat2cell([1:length(validC)]', ones(length(validC),1), 1);
C = cat(2, cellNum, validCellID);
T = cell2table(C, 'VariableNames', {'ID', 'Cell'});
writetable(T, sprintf('%sCells.txt', nameSubjNeural), 'Delimiter', '\t')