% genFigs_cellFeatureDataMining.m
%
% 2017/10/09 SHP
% Wondering around to find out something in cells' responses to movies
% in relation to the features or eye gaze (see also
% genFigs_eyeGazeDataMining.m)

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS';
    dirProcdata = '/Volumes/PROCDATA';
    dirDataHome = fullfile(dirProcdata, 'parksh');
    dirLibrary = '/Volumes/LIBRARY';
    %         addpath(fullfile(dirLibrary, 'matlab_utils'));
    dirDataNeural = fullfile(dirDataHome, 'Spi');
else % on virtual machine
    dirProjects = '/projects';
    dirProcdata = '/procdata';
    dirDataHome = fullfile(dirProcdata, 'parksh');
    dirLibrary = '/library';
    %         addpath(fullfile(dirLibrary, 'matlab_utils'));
    dirDataNeural = fullfile(dirDataHome, 'Spi');
end

%% Get the time series first

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
setMovie = [1 2 3];
% nt = 375;
% 
% % MION function
% TR=2.4;
% k = gampdf([-40:TR:40],4,2);

matSDF = struct([]);
curUnit = 0;
for iSubj = 1:length(setNameSubjNeural)
    
    % Load the parameter
    nameSubjNeural = setNameSubjNeural{iSubj};
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
    load(fullfile(dirDataNeural, filenameNeural), 'paramSDF')
    fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
        
    switch lower(nameSubjNeural)
        case 'spi'
            excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
            validC = setdiff(1:length(paramSDF.setCellIDs), excChanIndex)';
        otherwise
            [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
            validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
    end   
    matSDF(iSubj).setCellIDs = paramSDF.setCellIDs(validC);
    
    indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
    % 3. Fine temporal resolution
    clear FR_dTfine
    FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
        setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
    % concatenate across movies
%     matFR_SU=cell(135, );
    for iUnit = 1:size(FR_dTfine,1)
        curUnit = curUnit + 1;
        tempFR = cat(2, FR_dTfine(iUnit, :).matFR); 
        matFR_SU(curUnit, :) = tempFR; 
        
        %cat(1, FR_dTfine(iUnit, :).mnFR);
        %matFR_SU(:,iUnit) = tempFR;
    end
%     matSDF(iSubj).matFR_SU = matFR_SU;
    
end

% for iCell = 1:size(matFR_SU, 1)
% for iMov = 1:size(matFR_SU, 2)
% %     numTrial(iCell, iMov) = size(FR_dTfine(iCell, iMov).matFR{1}, 2);
% avgSTD(iCell, iMov) = mean(std(matFR_SU{iCell, iMov}'));
% end
% end

for iCell = 1:size(matFR_SU, 1)
for iMov = 1:size(matFR_SU, 2)
tempR = corrcoef(matFR_SU{iCell, iMov});
avgR(iCell, iMov) = mean(tempR(tril(tempR,-1)>0));
end
end

% 1) Clustering based on corr maps
load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  


%% 
setK = paramClustering_global.setK; %Clustering.setK;

matWSS_corrMap=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS_corrMap(:,iK) = sum(Clustering_moviemask_valid.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
end

totalSS_corrMap = Clustering_moviemask_valid.totalSS_SU;
% betweenSS_corrMap = totalSS_corrMap-matWSS_corrMap;
propExplained_corrMap = (totalSS_corrMap-matWSS_corrMap)./totalSS_corrMap; %matExpVar./totalSS;


%% Compare different clustering
K=7;
locMode_corrMap = find(propExplained_corrMap(:,K-1)==mode(propExplained_corrMap(:,K-1)));
indClust_SU = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode_corrMap(1)); % based on corr map
[sortedClust, indSortChan]=sort(indClust_SU);

oldIndCluster = [1 2 5 6 3 7 4]; % [3 4 5 2 1]; % [4 1 6 3 5 2 7]; 
indSortChan_new = [];
for iC = 1:K %7
    curC = oldIndCluster(iC);
    tempind = indSortChan(sortedClust==curC);
    indSortChan_new = cat(1, indSortChan_new, tempind);
end
sortedClust_new = indClust_SU(indSortChan_new);

xPos = find(abs(diff(sortedClust_new))>0)+0.5;
line([xPos xPos]', repmat([0; 37], 1, 6), 'Color', 'k') 
 
