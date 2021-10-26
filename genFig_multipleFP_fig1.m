% FIg1B
% example corr map
nameSubj = 'Tor';  %'Spi';
load(sprintf('/procdata/parksh/_macaque/%s/CorrMap_SU_%sArtMovie123_new.mat', nameSubj, nameSubj));

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';
nx = 40; ny = 64; nz = 32; %LR-AP-DV

indUnit = 7; %6; %17; %9;
matR_3d = reshape(matR_SU(:, indUnit), [nx, ny, nz]);

fig_corrMap = figure;
set(fig_corrMap, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [121   643   615   386])
iSlice = 7;
imagesc(squeeze(matR_3d(iSlice, :, :))')
colormap(gray)
title(sprintf('Cell %d: Slice %d', indUnit, iSlice))
axis off
set(gca, 'CLim', [-1 1].*.4)
print(fig_corrMap, fullfile(dirFig, sprintf('multipleFP_Fig1B_exampleCorrMap_%sCell%dSlice%', nameSubj, indUnit, iSlice)), '-r200', '-dtiff');

for iSlice = 1:nx
    figure(fig_corrMap);
    imagesc(squeeze(matR_3d(iSlice, :, :))')
    colormap(gray)
    title(sprintf('Cell %d: Slice %d', indUnit, iSlice))
    input('')
end

