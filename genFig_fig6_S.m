% genFig_fig6_S.m
%
% Correlation between the time series of Cell Group 3, High-gamma LFP (from three
% monkeys), and AF seed voxel (from two monkeys) that are used in Figure 6

%% Settings
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
dirDataHome = fullfile(dirProcdata, 'parksh');
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/'); % for saving figures as graphic files

% Add necessary toolbox  
addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution
addpath(fullfile(dirProjects, 'parksh/NeuralBOLD/analysis/BlockAna/BERscripts/')) % for 'decimate3D'

setNameSubjNeural = {'Tor', 'Rho', 'Sig'}; %'Sig'; %'Rho'; %'Tor';
setNameSubjBOLD = {'Art', 'Ava'}; %'Sig'; %'Rho'; %'Tor';


%% Load data and make the time course
% LFP
matLFP = zeros(375, length(setNameSubjNeural));

for iSubj=1:length(setNameSubjNeural)
    nameSubjNeural = setNameSubjNeural{iSubj};
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    filenameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
    
    load(fullfile(dirDataNeural, filenameNeural_BLP))
    
    
    % % 1. gamma pdf
    TR=2.4;
    k = gampdf([-40:TR:40],4,2);
    iF = 4; % high-gamma frequency
    setMovie = [1 2 3];
    neuralrgrs=[];
    for iMov = 1:length(setMovie)
        curNeuralTC = BLPRGR(setMovie(iMov)).meanBLP{iF}(8:125); % S(validC(iChan), indMovieNeuron(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
        curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
        curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
        curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
        
        neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
    end
    
    matLFP(:,iSubj) = neuralrgrs';
end


% Seed Voxel
matSeedVoxel = zeros(375, length(setNameSubjBOLD)*2); % For each subject, mean ROI and one voxel

for iSubj=1:length(setNameSubjBOLD)
    nameSubjBOLD = setNameSubjBOLD{iSubj};
    dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);
    filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];
    
    load(fullfile(dirDataBOLD, filenameBOLD))

    % 1. fMRI tc
    fmritc=[];
    % indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
    for iM = 1:length(setMovie)
        curvoltc = voltcIndMov{setMovie(iM)}; % flip the MION signal 
        avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
        if ~isempty(find(avgvoltc==0, 1))
            avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
        end
        pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
        fmritc = cat(4,fmritc,pcvoltc);
    end
    fmritc = fmritc.*(-1);
    
    % Get the ROI
    dirROI = fullfile(dirDataBOLD, 'ROIs');
    d_face = dir(fullfile(dirROI, '*faceROIs2.mat'));
    load(fullfile(dirROI, d_face.name));
    tempName_faceROIs = char(fieldnames(load(fullfile(dirROI, d_face.name))));
    eval(['faceROIs=', tempName_faceROIs, ';']) % get face ROI data into "faceROIs"
    eval(['clear ' tempName_faceROIs]) % clear up

    % Get ROI names & coordinates in EPI coords
    resizeFactor = [3 3 3]; % DSP.proc.params3d.res./ROIs(1,1).params.res; % calculate the scaling factor
    
    switch lower(nameSubjBOLD)
        case 'art'
            indAF = 4;
        case 'ava'
            indAF = 7;
    end
    nameROI = faceROIs(indAF).name(strfind(faceROIs(indAF).name, '_')+1:end); %faceROIs(iFR).name(strfind(faceROIs(iFR).name, '_')+1:end);
    voxROI=decimate3D(faceROIs(indAF).vol3D, resizeFactor, .25); %decimate3D(faceROIs(iFR).vol3D, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs
    [a, b, c] = ind2sub(size(voxROI), find(voxROI==1)); % Get indices of AF voxels in EPI 3D coords
    indVox_ROI = [a b c];
    indVox_ROI_sub = sub2ind(size(voxROI), a, b, c);
    
    matBOLD_ROI=[];
    for iVox = 1:length(a)
        matBOLD_ROI(:,iVox) = fmritc(a(iVox), b(iVox), c(iVox),:); %matBOLD_shuffle(a(iVox), b(iVox), c(iVox),:); %matBOLD(a(iVox), b(iVox), c(iVox),:);
    end
    meanBOLD_ROI = nanmean(matBOLD_ROI, 2);
    % steBOLD_ROI = nanstd(matBOLD_ROI, [], 2)./sqrt(size(matBOLD_ROI,2)-1);
    
    iVox=2; %8;
    curRGR = matBOLD_ROI(:,iVox);% meanBOLD_ROI;
    
    matSeedVoxel(:,(iSubj-1)*2+1) = meanBOLD_ROI;
    matSeedVoxel(:,(iSubj-1)*2+2) = curRGR;

end


% Cluster 3 from Toroid
nameSubjNeural = 'Tor';
nameSubjBOLD = 'Art';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);

load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)));
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr')
load(fullfile(dirDataNeural, sprintf('%s_movieTS_SU_indMov.mat', nameSubjNeural)))

matIndClust_SU = cat(2, Clustering.resultKMeans.SU_indCluster);
curK = 7; 
[sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));

setMovie = [1 2 3];
validC = paramCorr.validChanIndex; %find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]

% MION function
TR=2.4;
k = gampdf([-40:TR:40],4,2);

matFR_TR = zeros(375, length(validC));
for iChan = 1:length(validC) 

    % Modified 2016/04/05, 2016/04/27 by SHP
    neuralrgrs=[];
    for iMov = 1:length(setMovie)
        curNeuralTC = S(validC(iChan), setMovie(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
        curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
        curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
        curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
        
        neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
    end
    
    matFR_TR(:,iChan) = neuralrgrs';
end

% Average SDF in each cluster across cells
matFR_cluster3=[];
oldIndCluster =  [4 1 6 3 5 2 7]; %
iK=3;
indClust = oldIndCluster(iK);
matFR_cluster3 = matFR_TR(:,indSortChan(sortedClust==indClust));
matFR_cluster3 = cat(2, matFR_cluster3, nanmean(matFR_cluster3, 2)); % add average TS in cluster 3

% Average SDF across all the units
meanFR_TR = mean(matFR_TR, 2);

%% Compute correlation across different time series
% Matrix of time series (375 x 16) in the order of...
%   - 1:3   : high-gamma LFP from M3, M4, M5, MION convolved
%   - 4:5   : mean AF ROI time series from M1, one AF voxel time series from M1
%   - 6:7   : mean AF ROI time series from M2, one AF voxel time series from M2
%   - 8:14 : all the seven units (082a, 090a, 095a, 100a, 101a, 105a, 122b) in
%   Cluster 3, MION convolved
%   - 15    : averaged time series across seven units in Cluster 3
%   - 16    : averaged time series across ALL 48 units
matTS = cat(2, matLFP, matSeedVoxel, matFR_cluster3, meanFR_TR); 

[matR_TS] = corr(matTS, 'rows', 'complete', 'type', 'Spearman');


%% Plot the results
fig6_S = figure;
set(fig6_S, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 800 800])

imagesc(tril(matR_TS))
axis off
rgb=hslcolormap(256, 'bc.yr', 1, [0.2 1 0.2]); colormap(rgb)
caxis([-1 1])

% save
print(fig6_S, fullfile(dirFig, 'fig6_S'), '-r150', '-dtiff')

% colorbar
figure(fig6_S);
colorbar;
print(fig6_S, fullfile(dirFig, 'fig6_S_colorbar'), '-r150', '-dtiff')
