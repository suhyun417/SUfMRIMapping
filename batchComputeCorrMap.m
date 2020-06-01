% batchComputeCorrMap.m
%
% 2020/05/12 SHP

clear all;

setNameSubjNeural = {'Tor', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dex'}; %, 'Rho', 'Dav'}
setNameSubjBOLD = {'Art', 'Ava'};

setMovie = [1 2 3];
flagSaveFile = 1;
typeMION = 1;

for iSN = 1:length(setNameSubjNeural)
    nameSubjNeural = setNameSubjNeural{iSN};
    
    for iSB = 1:length(setNameSubjBOLD)
        nameSubjBOLD = setNameSubjBOLD{iSB};
        
        % compute correlation between each neuron and all the voxels
        computeCorrMap(nameSubjNeural, nameSubjBOLD, setMovie, flagSaveFile, typeMION);
        
    end
end

