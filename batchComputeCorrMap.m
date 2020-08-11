% batchComputeCorrMap.m
%
% 2020/05/12 SHP
% 2020/07/18 SHP: compute pca-res maps


clear all;

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dav'}; %'Dex'}; %, 'Dav'}
setNameSubjBOLD = {'Art', 'Ava'};

setMovie = [1 2 3];
flagSaveFile = 1;
typeMION = 1;

for iSB = 1:length(setNameSubjBOLD)
nameSubjBOLD = setNameSubjBOLD{iSB}; %'Art'; %setNameSubjBOLD{iSB};

computeCorrMap_pcares_masked('Dav', nameSubjBOLD, setMovie, flagSaveFile);
for iSN = 1:length(setNameSubjNeural)
    nameSubjNeural = setNameSubjNeural{iSN};
    
    %         % compute correlation between each neuron and all the voxels
    %         computeCorrMap(nameSubjNeural, nameSubjBOLD, setMovie, flagSaveFile, typeMION);
    
    % compute correlation between each neuron and all the voxels: pcares
    fprintf(1, '      ::     Computing pca-res corr map for %s (%d/%d) x %s  (%s)     :: \n\n ', nameSubjNeural, iSN, length(setNameSubjNeural), nameSubjBOLD, datestr(now, 'mm/dd/yy HH:MM:SS'))
    computeCorrMap_pcares(nameSubjNeural, nameSubjBOLD, setMovie, flagSaveFile);
        
end
end


% nameSubjBOLD = 'Ava'; %setNameSubjBOLD{iSB};
% 
% for iSN = 1:length(setNameSubjNeural)
%     nameSubjNeural = setNameSubjNeural{iSN};
%     
%     %         % compute correlation between each neuron and all the voxels
%     %         computeCorrMap(nameSubjNeural, nameSubjBOLD, setMovie, flagSaveFile, typeMION);
%     
%     % compute correlation between each neuron and all the voxels: pcares
%     computeCorrMap_pcares(nameSubjNeural, nameSubjBOLD, setMovie, flagSaveFile);
%     
%     %         if sum(strcmpi(nameSubjNeural, {'dan', 'dav'}))>0
%     computeCorrMap_pcares_masked(nameSubjNeural, nameSubjBOLD, setMovie, flagSaveFile);
%     %         end
%     
% end
