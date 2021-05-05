% genFig_multipleFP_S_waveform.m
%
% Generate a supplementary figure for multiple face patch paper,
% which shows individual single unit correlation map
% for two example neurons for units from the same channel
% Draw waveforms for Fig 2c units: modified from genFig_fig2.m and
% genFigs_waveformDataMining.m
%
% 2021/02/23 SHP
% 2021/05/03 SHP: modified the FSI computation part


% Directory for saving figures as graphic files
dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs/';


%% Which units?
setExampleCellIDs = {'29Dav', '30Dav', '130AFMoc', '131AFMoc', '28AMWas', '29AMWas', '108AMMoc', '109AMMoc'};

%%
figSWaveForm = figure; % for each unit
set(figSWaveForm, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 100 300])

for iCell = 1:length(setExampleCellIDs)
    curCellID = setExampleCellIDs{iCell};
    nameSubjNeural = char(curCellID(end-2:end));
    
%     dirDataNeural = fullfile(dirProcdata, '/parksh/_macaque', nameSubject);
%     
%     setMovIDs = [1 2 3];
    
    d_n = dir(sprintf('/procdata/parksh/_macaque/%s/*sig*.mat', nameSubjNeural));
    fname = {d_n.name}';
    
    cellID = char(curCellID(1:end-3));
    tLoc = find(contains(fname, cellID)>0, 1);
    load(fullfile(d_n(tLoc).folder, d_n(tLoc).name))
    catWF = cat(2, dat.waveform{:});
    
    tAxis_ms = (1:size(catWF,1)).*(1/24414.0625)*1000;
    
%     figure;
%     plot(tAxis_ms, mean(catWF,2), 'k-', 'LineWidth', 3)
%     title(sprintf('Cell %s:', curCellID))
    
    figure(figSWaveForm); clf;
    subplot('Position', [0 0 1 1]);
    plot(tAxis_ms, mean(catWF,2), 'k-', 'LineWidth', 3)
    ylim([-80 25])
    xlim([0.5 2.5])
    axis off

    % save
    saveFileName = sprintf('multipleFP_FigS_waveform_%02d_%s', iCell, curCellID);
%     print(figSWaveForm, fullfile(dirFig, saveFileName), '-depsc')
    print(figSWaveForm, fullfile(dirFig, saveFileName), '-r200', '-dtiff')

end

%% Scale bar
bar_time_ms = 1; % scale bar for time: 1 ms
bar_amp_microV = 20; % scale bar for amplitude: 20 micro volt

figSWaveForm_scale = figure; 
set(figSWaveForm_scale, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 100 300])
clf
subplot('Position', [0 0 1 1]);
line([1 1; 1 1+bar_time_ms], [-50 -50; -50+bar_amp_microV -50], 'Color', 'k', 'LineWidth', 3)
ylim([-80 25])
xlim([0.5 2.5])
axis off
 % save
 saveFileName_sb = 'multipleFP_FigS_waveform_scaleBar_1ms_20uv';
%  print(figSWaveForm_scale, fullfile(dirFig, saveFileName_sb), '-depsc')
 print(figSWaveForm_scale, fullfile(dirFig, saveFileName_sb), '-r200', '-dtiff')
 
 
 %% Face selectivity of these examples
 setExampleCellIDs = {'29Dav', '30Dav', '130AFMoc', '131AFMoc', '28AMWas', '29AMWas', '108AMMoc', '109AMMoc'};

 % Read the spreadsheet to load corresponding fingerprinting data file name
filename_xls = '/procdata/parksh/_macaque/multipleFP_4FPneurons_CellIDFingerPrinting.xls';
C = readcell(filename_xls); % 1st col: cell ID in movie data, 2st col: fingerprinting results directory, 3rd col: cell file name

%% Cell-by-cell computation of FSI etc.
cond_face = 1:2; % human face & monkey face
cond_obj = 4; % object

matFaceSelective = NaN(size(setExampleCellIDs));
for iCell = 1:numel(setExampleCellIDs)
    curCellID = setExampleCellIDs{iCell};
    nameSubjNeural = char(curCellID(end-2:end));
    
    % get the index from the spreadsheet 
    indCell = find(strcmp(C(:,1), curCellID)>0);
    
    if ~cellfun(@ischar, C(indCell,2)) % has zero (not character) if there's no fingerprinting results        
        matFaceSelective(iCell) = 0.5; % NaN for gray in the furture plot
        continue;
    end
    
    % load fingerprinting file
    load(sprintf('/procdata/parksh/_macaque/%s%s/%s.mat', nameSubjNeural, C{indCell, 2}, C{indCell, 3}))
    
    % compute fsi & other things
    fprintf(1, ':::::%s:::::\n', curCellID)
    curFSI = calc_fr_from_multidays_for_fsi(multiday, cond_face, cond_obj);
    
    if curFSI > 0.33 % abs(curFSI) > 0.33
        matFaceSelective(iCell) = 1; %face-selective
    else
        matFaceSelective(iCell) = 0; % not face-selective
    end

end


%% Quick overview of Fig 2b
figure;
imagesc(matFaceSelective);
colormap(gray)
hold on;

% T = text(xC(:), yC(:), setExampleCellIDs(:), 'HorizontalAlignment', 'center');
set(T(matFaceSelective<0.5), 'Color', 'w')
set(gca, 'XTick', [], 'YTick', []);

