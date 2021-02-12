% genFig_multipleFP_fig2_FSI.m
%
% 2021/02/12 SHP
% Compute face selectivity for example neurons shown in Figure 2

%% 2021/02/06 NOTES from "computeFSI.m"
% Start from the spreadsheet to go back to the original cell ID of channel
% number + alphabet coding
% Once you get the original cell ID, load the FPrint results (.mat file) of 
% that cell from FPrint directory
% Compute firing rate for each stimulus, normalize them, then compute face selectivity index
%       -- refer to "response_calc_SHP" for basic concatenation of multiday
%       data
%       -- refer to "plot_heatmap_FPrint" for normalization & FSI
%       computation
% Put the results back in the matrix of neurons x stimulus, ideally
% neurons in the order of correlation map order of 389 cells (refer to 
% "computeCorrMap_AllCells_FPAreaSummary") to make it easier to link them
% to the correlation map results

%% Figure 2 example cells
setExampleCellIDs = {'33Dav', '130AFMoc', '097aMat', '10Dan', ...
    '27Dav', '065aTor', '39AMWas', '117AMMoc', ...
    '25Dav', '022bSpi', '51AMWas', '109AMMoc', ...
    '16Dav', '045aSpi', '33AMWas', '05Dan', ...
    '06Dav', '122AFMoc', '06AMWas', '115AMMoc'};

