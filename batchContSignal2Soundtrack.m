% batchContSignal2Soundtrack.m
Script;

% for Principal components
dirData = '/procdata/parksh'; %'/Volumes/PROCDATA/parksh';
dirDataNeural = '/procdata/parksh/Tor/';
dirDataEig = '/procdata/parksh/Tor/eigen/';
dirDataEigSound = [dirDataEig, 'audioPC'];

setMovID = [1 2 3]; %[1 2 3 10 11 12 13 14 15];
validC = [1:8, 11:34, 36:51]; % valid channel with movie [1 2 3]

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

for iPC = setPC
    
    saveFileName = fullfile(dirDataEigSound, sprintf('PC%d.wav', iPC));
    inputSig = pcaMovCat_1sec.coeff(:,iPC); %coeff(:,iPC);
    Fs_org = 1;
    flagAM = 1; %amplitude modulation
    contSignal2Soundtrack(inputSig, Fs_org, flagAM, saveFileName)
    
end