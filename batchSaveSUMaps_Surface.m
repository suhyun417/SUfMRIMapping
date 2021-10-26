% batchSaveSUMaps_Surface.m
%
% 2020/07/21 SHP


clear all;

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dav'}; %, 'Dex'}; %, 'Dav'}
setNameSubjBOLD = {'Art', 'Ava'};

setMovie = [1 2 3];

for iSB = 1:length(setNameSubjBOLD)
    nameSubjBOLD = setNameSubjBOLD{iSB};
    
    for iSN = 1:length(setNameSubjNeural)
        nameSubjNeural = setNameSubjNeural{iSN};
        
        % Saving SU maps (unmasked pca-res maps)
        fprintf(1, '      ::     Saving pca-res corr map to surface for %s (%d/%d) x %s  (%s)     :: \n\n ', nameSubjNeural, iSN, length(setNameSubjNeural), nameSubjBOLD, datestr(now, 'mm/dd/yy HH:MM:SS'))
        saveSUMaps_Surface_pcares(nameSubjNeural, nameSubjBOLD, setMovie);
        
        cd /projects/parksh/NeuroMRI/analysis
    end
end
