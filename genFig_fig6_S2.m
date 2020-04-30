% genFig_fig6_S2.m
%
% LFPs in different frequency bands and their consistency across subjects
% Another supplementary for Figure 6

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

% % Add necessary toolbox  
% addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution
% addpath(fullfile(dirProjects, 'parksh/NeuralBOLD/analysis/BlockAna/BERscripts/')) % for 'decimate3D'

setNameSubjNeural = {'Tor', 'Rho', 'Sig'}; %'Sig'; %'Rho'; %'Tor';
% setNameSubjBOLD = {'Art', 'Ava'}; %'Sig'; %'Rho'; %'Tor';


%% Load data and make the time course
% LFP
matLFP = zeros(375, length(setNameSubjNeural), 4);

for iSubj=1:length(setNameSubjNeural)
    nameSubjNeural = setNameSubjNeural{iSubj};
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    filenameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
    
    load(fullfile(dirDataNeural, filenameNeural_BLP))
    
    cc = cat(1, BLPRGR(1:3).matBLP);
    
    figure(100+iSubj);
    for iFR=1:4
        for iM=1:3
        subplot(4,1,iFR)
        plot(125*(iM-1)+1:125*iM, cc{iM, iFR}); hold on;
        end
    end
    
    tempC = cat(1, BLPRGR(1:3).meanBLP);
    tempM = cell2mat(tempC);
    
    matLFP(:,iSubj,:) = tempM;
end
    
for iSubj=1:length(setNameSubjNeural)
    nameSubjNeural = setNameSubjNeural{iSubj};
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    filenameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
    
    load(fullfile(dirDataNeural, filenameNeural_BLP))
    % % 1. gamma pdf
    TR=2.4;
    k = gampdf([-40:TR:40],4,2);
    for iF =1:4; % high-gamma frequency
    setMovie = [1 2 3];
    neuralrgrs=[];
    for iMov = 1:length(setMovie)
        curNeuralTC = BLPRGR(setMovie(iMov)).meanBLP{iF}(8:125); % S(validC(iChan), indMovieNeuron(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
        curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
        curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
        curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
        
        neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
    end
    
    matLFP_MION(:,iSubj,iF) = neuralrgrs';
    end
end



%% Plot the results
fig6_S2 = figure;
set(fig6_S2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 800 800])

imagesc(tril(matR_TS))
axis off
rgb=hslcolormap(256, 'bc.yr', 1, [0.2 1 0.2]); colormap(rgb)
caxis([-1 1])

% save
print(fig6_S2, fullfile(dirFig, 'fig6_S2'), '-r150', '-dtiff')

% % colorbar
% figure(fig6_S2);
% colorbar;
% print(fig6_S2, fullfile(dirFig, 'fig6_S2_colorbar'), '-r150', '-dtiff')
