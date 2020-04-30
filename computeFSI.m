% computeFSI.m
% 
% 2019/03 SHP
%       - Compute face-selective index using responses to fingerprinting stimulus set
%       for neurons that Kenji Koyano and Elena Waidmann
%       - Modified mainly from "/projects/parksh/_toolbox/plot_heatmap_FPrint.m" written originally by Kenji Koyano
%       - Cell selection part is copied from "./renameFiles_Spice.m"


%%%%%% Still need to put the pieces together %%%%%%%%%%%%

%% Cell Selection (read the xls file, choose the cells that have positive "evaluation" values)
clear all;

nameSubjNeural = 'Spi';
directory.dataHome = '/procdata/parksh';

% Directory info
directory.dataNeural = fullfile(directory.dataHome, nameSubjNeural); % '/procdata/parksh/Spi';
% dirDataNeural_org = fullfile(dirDataNeural, '_orgData'); % '/procdata/parksh/Spi/_orgData';
directory.source = fullfile(directory.dataNeural, '_orgData', 'Spice180120_26_movie');
directory.destination = fullfile(directory.dataNeural,  '2018Jan_movie');

% nameSession_source = 'Spice180120_26_movie';
% nameSession_destination = '2018Jan_movie';

% read the excel file of cell list
filename_xls =  fullfile(directory.dataNeural, '_orgData/Spice180120_other/Spice_cell_list_180120_SHP.xls');
T = readtable(filename_xls, 'ReadRowNames', true);

 if T{indCell,5} > 0 % if the evaluation is okay
 end
 
 %% From "genFigs_multiplePatchDataMining.m"
 %% Fingerprinting results
 dirDataHome = '/procdata/parksh';
setNameSubjNeural = {'Ava', 'Dav', 'Spi', 'Mat', 'Dan'}; % {'Tor', 'Rho', 'Sig', 'Spi'};
setMovie = [1 2 3]; %[4 5 6]; %[1 2 3];
tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));

iSubj = 5; %5; %2; %5;

nameSubjNeural = setNameSubjNeural{iSubj};
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
switch lower(nameSubjNeural)
    case 'spi'
        nameSession_FPrint = 'Spice180120_other';
    case 'dav'
        nameSession_FPrint = 'Davida180515_other';
    case 'dan'
        nameSession_FPrint = 'Dango180123_other';
end
dirFPrint = fullfile(dirDataNeural, '_orgData', nameSession_FPrint, 'FPrint');

% get the file names of finger printing data
clear filename_fp cellID_fp d_fp
d_fp = dir(fullfile(dirFPrint, '*.mat'));
[filename_fp{1:length(d_fp)}] = deal(d_fp.name);
for iF = 1:length(filename_fp)
cellID_fp{iF} = filename_fp{iF}(1:strfind(filename_fp{iF}, '_')-1);
end

% get the list of selected cells
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%sArtMovie%s_new.mat', nameSubjNeural, MovieStr)), 'paramCorr')

% going through relevant cells and get the finger printing results
matFR_fp = NaN(10, 6, size(paramCorr.validChanID, 1));
clear catDate
for iCell = 1:size(paramCorr.validChanID, 1)
    tt = strcmp(sprintf('%d', str2double(paramCorr.validChanID(iCell,:))), cellID_fp);
    indCell = find(tt);
    [curMat nDate] = response_calc_SHP(fullfile(dirFPrint, filename_fp{indCell}));
    matFR_fp(:, :, iCell) = curMat;
    catDate(iCell, 1) = nDate;
end
clear catMat
for iCell = 1:size(matFR_fp, 3)
    curMat = matFR_fp(:,:,iCell);
    curMat_norm = curMat./max(max(curMat));
    catMat(:,iCell) = curMat_norm(:);
end
figure
imagesc(catMat')
colormap(hot);
colorbar;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
title(sprintf('%s neurons: Fingerprinting normalized responses', nameSubjNeural))
ylabel('Cells')
xlabel('Stimulus')
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
print(gcf, fullfile(dirFig, sprintf('%s_fingerprinting_Movie%s', nameSubjNeural, MovieStr)), '-r150', '-dtiff')

matFR_fp_mean = squeeze(mean(matFR_fp));
figure; bar(matFR_fp_mean')

matFR_fp_faceother=[];
matFR_fp_faceother(1,:) = mean(matFR_fp_mean(1:2, :));
matFR_fp_faceother(2,:) = mean(matFR_fp_mean(3:6, :));
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1400 300])
bar(matFR_fp_faceother')
set(gca, 'XTick', 1:length(matFR_fp_faceother))
xlabel('Cell Index')
title(sprintf('%s: fingerprinting results: movie %s: faces vs. others', nameSubjNeural, MovieStr));
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
print(gcf, fullfile(dirFig, sprintf('%s_fingerprinting_Movie%s_avgFacesOthers', nameSubjNeural, MovieStr)), '-depsc')

 
%% FSI computation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% needs to be modified for "multiday" merged FPrint results

% /procdata/parksh/Spi/_orgData/Spice180120_other/FPrint
% /procdata/parksh/Dav/_orgData/Davida180515_other/FPrint
% /procdata/parksh/Dan/_orgData/Dango180123_other/FPrint

% multiday = 
% 
%        response: [5x60 double]
%            date: {'180515'  '180516'  '180517'  '180518'}
%      rasterdata: {[651x643 double]  [651x734 double]  [651x669 double]  [651x665 double]}
%          offset: {[643x1 double]  [734x1 double]  [669x1 double]  [665x1 double]}
%            stim: {[643x1 double]  [734x1 double]  [669x1 double]  [665x1 double]}
%        stimlist: [60x1 double]
%              wf: [62x4 double]
%            wfsd: [98x4 double]
%       rasterwin: [-200 450]
%        dateindx: [4x4 logical]
%              fr: [1x1 struct]
%            corr: [1x1 struct]
%          fr_win: [70 350]
%     fr_base_win: [-100 50]

TankName = 'Mochi181128';
directory.session = ['/procdata/koyanok/physiology/raster/' TankName '/FPrint1/'];
directory.save    = ['/procdata/koyanok/physiology/raster/' TankName '/FPrint1/'];
% directory.session = ['/procdata/waidmannen/physiology/raster/' TankName '/FPrint1/'];
% directory.save    = ['/procdata/waidmannen/physiology/raster/' TankName '/FPrint1/'];
% directory.save    = ['/procdata/koyanok/physiology/analysis/sessions/Matcha180621/heatmaps/'];
response_win = [50 150];   % time window for response firing rate calculation, ms
baseline_win = [-100 50];  % time window for baseline firing rate calculation, ms
plot_order = [1:60];  % order of stimuli for plotting
category_boundary = 0:10:60;    % boundary of the category
category_txt = {'HFace';'MFace';'MBody';'Object';'Scene';'Bird'}; % label for stimulus category, after reordering

%% list cells
celllist = dir([directory.session '*.mat']);
celllist = {celllist.name};
% tmpcelllist = dir([directory.session '*.mat']);
% tmpcelllist = {tmpcelllist.name};
% celllist = cell(1,61);
% for i=41:101
% 	celllist{i-40} = tmpcelllist{i};
% end
%% loop for cell
norm_fr =[];    % normalized firing rate
for cellno = 1:length(celllist)
    [~,cellname,~] = fileparts(celllist{cellno});
    disp(cellname);
    
    cellnamelist{1,cellno} = cellname;
    
    %% load data and format raster data into an array
    load([directory.session cellname '.mat']);
    
    %% calculate firing rate of each trial
    response_calculate_win = response_win - raster.win(1);
    fr_trial = sum(raster.data(response_calculate_win(1):response_calculate_win(2),:),1)/diff(response_win)*1000;  % Hz
        
    %% loop for stimulus
    fr = nan(length(code.event),1);
    for stimno = 1:length(code.event)
        stimcode = code.event(stimno);
        stim_indx = raster.stimcode==stimcode;
        fr(stimno) = nanmean(fr_trial(stim_indx));
    end
    fr = fr(plot_order);    % sort stimuli
        
    %% normalize firing rate
    baseline_calculate_win = baseline_win - raster.win(1);
    baseline_fr_trial = sum(raster.data(baseline_calculate_win(1):baseline_calculate_win(2),:),1)/diff(baseline_win)*1000;  % Hz
    baseline_fr = nanmean(baseline_fr_trial);
%     max_fr = nanmax(fr);
    max_fr = nanmax(abs(fr-baseline_fr));
    norm_fr(:,cellno) = (fr - baseline_fr)/(max_fr);
    
    %% calculate face-selective index
    fsi_fun = @(x,y) (x-y)/(abs(x)+abs(y))*sign(double(sign(x)>0|sign(y)>0)-0.5);
    response_face    = nanmean(fr(1:30))-baseline_fr;    % Human face, monkey face and monkey body
%     response_nonface = nanmean(fr(40:60))-baseline_fr;   % Object, scene and bird
%     response_face    = nanmean(fr(1:20))-baseline_fr;    % Human face and monkey face
    response_nonface = nanmean(fr(40:50))-baseline_fr;   % Object and scene
    fsi(cellno,1) = fsi_fun(response_face,response_nonface);
end

%% save data
fsidata.fsi = fsi;
fsidata.norm_fr = norm_fr;
fsidata.info.cellname = cellnamelist;
fsidata.info.duration.response = response_win;
fsidata.info.duration.baseline = baseline_win;

if ~exist(directory.save)
    mkdir(directory.save);
end
save([directory.save 'fsidata_AF' TankName '.mat'],'fsidata');


%% sort normalized firing rate matrix
sort_indx = [fsi,(1:length(fsi))'];
sort_indx = sortrows(sort_indx);
sort_indx = sort_indx(:,2);

%% plot population heat map, sorted by fsi
figure('Position',[200 200 600 400],'Color','w');
hold on;
imagesc(norm_fr(:,sort_indx)');
for ii = 1:length(category_boundary)
    plot([category_boundary(ii),category_boundary(ii)]+0.5,[0, size(norm_fr,2)]+0.5,'w');
end
set(gca,'XTick',(category_boundary(1:end-1)+diff(category_boundary)/2)+0.5);
set(gca,'XTickLabel',category_txt);
set(gca,'YTick',1:length(celllist),'YTickLabel',cellnamelist(sort_indx));
xlabel('Stimuli');
ylabel('Cells');
xlim([0 size(norm_fr,1)]+0.5);
ylim([0 size(norm_fr,2)]+0.5);
title([TankName ', Finger Printing (n=' num2str(length(celllist)) ')'],'Interpreter','none','FontSize',14);
colormap('hot');
colorbar;

print(gcf,[directory.save 'heatmap_population' TankName '.eps'],'-depsc2');
     
     