% genFigs_pcaDataMining


%% Principal Components and Movie contents
dirData = '/procdata/parksh'; %'/Volumes/PROCDATA/parksh';
dirDataNeural = '/procdata/parksh/Tor/';
dirDataEig = '/procdata/parksh/Tor/eigen/';

setMovID = [1 2 3]; %[1 2 3 10 11 12 13 14 15];
validC = [1:8, 11:34, 36:51]; % valid channel with movie [1 2 3]


% Use 1kHz SDF, Gaussian convolved
load(fullfile(dirData, 'dataNeuralBOLD_TorArt_indMov.mat'), 'S') % load 1kHz SDF

matSDF=[];
for iMov = 1:length(setMovID)
    movID = setMovID(iMov);
    tempMat=[];
    tempMat = cat(2,S(validC,movID).mnsdf);
    matSDF = [matSDF ; tempMat]; % concatenate across movies
end
matSDF_norm = matSDF-repmat(mean(matSDF), size(matSDF, 1), 1);
matSDF_norm = matSDF_norm./repmat(std(matSDF), size(matSDF,1), 1);
[coeff,pcascore,latent,tsquared,explained] = pca(matSDF_norm'); % coeff: time x PC (900000x47)

setPC = 1:5;

figure
clf
subplot(2,1,1)
imagesc(coeff(1:1000:900000, setPC)')
subplot(2,1,2)
for iP=setPC
    plot(coeff(1:1000:900000,iP)-(0.005*(iP-1)), '-')
    hold on
end
set(gca, 'YTick', -0.02:0.005:0, 'YTickLabel', 5:1)
set(gca, 'XMinorTick', 'on')


% Find when the PC was high (i.e. exceeds certain criterion in z-score)
resmpCoeff = coeff(1:1000:900000, setPC);
zscoreCoeff = zscore(resmpCoeff);

critHigh =1.5; %2; % in z-score
locHigh = find(zscoreCoeff>critHigh);
[i,j] = ind2sub(size(zscoreCoeff), locHigh);
taxis = 1:900; % in second
for iP=setPC
    loc =[];
    loc = i(j==iP);
    if exist('duration')
        timeHighPC{iP} = duration(0,0, taxis(loc)); % in minutes & seconds, to make it easy to compare with movies
    else % some versions of matlab doesn't have this function
        timeHighPC{iP} = taxis(loc); % in seconds
    end
    
end

locHigh_logical = (zscoreCoeff>critHigh);

% plot timepoints that exceed threhold
figure;
imagesc(zscoreCoeff(locHigh_logical)')
imagesc(locHigh_logical'.*(-1)+1)


% Use TR-resolution firing rate
load(fullfile(dirDataEig, 'pcaCatMovie123_FR_tor.mat'))
setPC=1:5;
zscoreCoeff_FR = zscore(pcaMovCat.coeff(:,setPC));
critHigh=1.5;
locHigh = find(zscoreCoeff_FR>critHigh);
[i,j] = ind2sub(size(zscoreCoeff_FR), locHigh);
axis_TR = 1.2:2.4:900;
for iP=setPC
    loc = [];
    loc = i(j==iP);
    if exist('duration')
        timeHighPC_FR{iP} = duration(0,0, axis_TR(loc)); % in minutes & seconds, to make it easy to compare with movies
    else % some versions of matlab doesn't have this function
        timeHighPC_FR{iP} = axis_TR(loc); % in seconds
    end    
end
locHigh_logical_FR = (zscoreCoeff_FR>critHigh);
imagesc(axis_TR, 1:5,locHigh_logical_FR'.*(-1)+1)
set(gca, 'YTick', 1:5, 'YTickLabel', 1:5)



%% Run PCA on 1sec data
d_n = dir(fullfile(dirDataNeural, '*sig*.mat'));
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files
% listMatSUFile = regexp(cat(2, d_n.name), '\w*.mat', 'match')'; %regexp(cat(2, d_n.name), '(?<=mat)\w*', 'match')';
% listSUchannelID = unique(regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match'));
listMovSU = unique(regexp(cat(2,d_n.name), '\d*(?=sig)', 'match'));
listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
listSUchannelID = unique(listSU_all);
setCellIDs = listSUchannelID; % should be cell array of string
setMovIDs = sort(str2num(char(listMovSU))');
S_1sec = createCellRegressor_indMov_discreteTime(dirDataNeural, setCellIDs, setMovIDs, 1);

matSDF_1sec=[];
for iMov = 1:length(setMovID)
    
    movID = setMovID(iMov);
    tempMat=[];
    
    tempMat = cat(2,S_1sec(validC,movID).mnFR); %cat(2,S(validC,movID).mnFR); %cat(2,S(validC,movID).mnsdf);
    matSDF_1sec = [matSDF_1sec ; tempMat]; % concatenate across movies
end


matSDF_norm_1sec = matSDF_1sec-repmat(mean(matSDF_1sec), size(matSDF_1sec, 1), 1);
matSDF_norm_1sec = matSDF_norm_1sec./repmat(std(matSDF_1sec), size(matSDF_1sec,1), 1);

[pcaMovCat_1sec.coeff, pcaMovCat_1sec.pcascore, pcaMovCat_1sec.latent, pcaMovCat_1sec.tsquared, pcaMovCat_1sec.explained] = pca(matSDF_norm_1sec');


%% Make a sound file out of PC
filename = 'pc1.wav';
workingDir_audio = '/Users/parks20/ResearchProjects/0NeuralBOLDCorr/movies/audioPC';
Fs=8000;
taxis = 1/Fs:1/Fs:900;
size(taxis)
resmpCoeff_1sec = resample(pcaMovCat_1sec.coeff(:,1), 8000, 1);
audiowrite(filename, resmpCoeff_1sec, Fs);
audioinfo(filename)


%
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

load('/procdata/parksh/Tor/eigen/pcaCatMovie123_SDF_tor.mat') %pcaCatMovie123_FR_tor.mat') %pcaCatMovie123.mat')

fNameHead = [dirFig 'pca_catMovie_'];

% PC time course
figPC = figure;
set(figPC, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 700 920 175])

for i=1:5
    figure(figPC); clf;
    plot(pcaMovCat.coeff(:,i), 'k-', 'LineWidth', 4)
    
    set(gca, 'LineWidth', 2, 'YTick', [], 'TickDir', 'out', 'XTickLabel', [])
    box off
    
    fName = sprintf([fNameHead 'timecourse_PC%d'], i);
    print(figPC, fName, '-depsc')
end

% PCA scores for each cell in 3d for clustering
modelCell = {'065a', '082a', '075a', '089b', '118a'};
indModelCell=[];
for i=1:length(modelCell)
    indModelCell(i) = find(strcmp(modelCell(i), paramPcaMov.validCellID));
end

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
title('PCA scores for each cell for top 3 PCs')
plot3(pcaMovCat.pcascore(:,1),pcaMovCat.pcascore(:,2),pcaMovCat.pcascore(:,3), 'ko',...
    'MarkerSize', 8, 'LineWidth', 2);
hold on;
plot3(pcaMovCat.pcascore(indModelCell,1),pcaMovCat.pcascore(indModelCell,2),pcaMovCat.pcascore(indModelCell,3), 'bo',...
    'MarkerSize', 8, 'LineWidth', 2)
grid on;
xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
text(pcaMovCat.pcascore(indModelCell,1),pcaMovCat.pcascore(indModelCell,2),pcaMovCat.pcascore(indModelCell,3),...
    paramPcaMov.validCellID(indModelCell,:), 'FontSize', 15, 'Color', 'b');

% PCA scores for each cell in colored plot
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
title('Coefficients of each PC')
imagesc(pcaMovCat.pcascore(:,1:5)) % only for the top 5 PCs
colorbar;
set(gca, 'YTick', 1:1:length(paramPcaMov.validCellID), 'YTickLabel', paramPcaMov.validCellID)
set(gca, 'XTick', 1:5)
xlabel('Principal component')
ylabel('Cell IDs')

% Percent variance explained
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
title('Explained Variance')
plot(cumsum(pcaMovCat.expVar), 'o', 'MarkerSize', 8, 'LineWidth', 3)
axis square
xlim([0 38]);
ylim([0 100])
set(gca, 'XTick', [10 20 30])
set(gca, 'LineWidth', 2, 'FontSize', 20)
box off
print(gcf, [fNameHead '_expVar'], '-depsc')

% first PC and face scale regressor (in TR unit)
figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(zscore(pcaMovCat.coeff(:,1)), 'b-')
% plot(zscore(resample(pcaMovCat.coeff(:,1), 0.1*10, 2.4*10)), 'bo-')
hold on
plot(scaleRGR, 'r-')
axis tight
legend('PC 1 (normalized)', 'Scale')
print(gcf, [fNameHead '_PC1andFaceScale'], '-depsc')

% Concatenated time course of each cell
matSDF=[];
for iMov = 1:length(pcaMovCat.movID)
    
    movID = pcaMovCat.movID(iMov);
    tempMat=[];
    
    tempMat = cat(2,S(paramPcaMov.validCell,movID).mnsdf);
    matSDF = [matSDF ; tempMat]; % concatenate across movies
end
% matSDF_norm = matSDF-repmat(mean(matSDF), size(matSDF, 1), 1);
% matSDF_norm = matSDF_norm./repmat(std(matSDF), size(matSDF,1), 1);

catSDF = matSDF(1:100:size(matSDF,1),:); % downsample
clear matSDF
figCatSDF = figure;
for iC = 1:size(catSDF,2)
    figure(figCatSDF); gcf;
    plot(catSDF(:,iC), '.')
    title(paramPcaMov.validCellID(iC,:))
    input('')
end




