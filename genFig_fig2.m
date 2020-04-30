% genFig_fig2.m
%
% Generate figure 2, which shows individual single unit correlation map
% for two example neurons (Fig 2a), for entire population (Fig 2b), and for
% units from the same channel (Fig 2c)
% Draw waveforms for Fig 2c units: modified from genFig_waveform.m
%
% 2016/04/19 SHP

%% Set the directory
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

nameSubject = 'Tor'; %'Rho'; % 'Tor'; %'Rho'; %'Tor';
dirNeuralData = fullfile(dirProcdata, '/parksh', nameSubject);
dirData_waveform = fullfile(dirNeuralData, '/_sortedSpikes/_waveformCompare');
dirData_spike = fullfile(dirNeuralData, '_sortedSpikes');
 

% toolkit for reading .plx file
addpath(fullfile(dirProjects, '/parksh/_toolbox/Waveform_Consistency-V1.1/'))

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');

%% Which units?
listChan = [82 85 86 89 91 96 122]'; 
unitLabel = {'a', 'b'};
numChan = length(listChan);

%% Read files and draw the waveform (for all channels)
fig2c = figure; % for each unit
set(fig2c, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 100 300])
for iChan = 1:numChan
    plxFileName = {};
    idChan = listChan(iChan);
    filename_wf = fullfile(dirData_waveform, sprintf('waveformComp_%d', idChan));
    load(filename_wf)
    
    % Option 1: Take a representative session
    curWF=[];
    curWF = data.waveform{1}; % some neurons don't have 2 units for the first session (don't know why..)
    if size(data.waveform{1},2)<2
        curWF = data.waveform{2};
    end
    
    % Plot
    
    for iUnit = 1:2
        figure(fig2c); clf;
        subplot('Position', [0 0 1 1]);
        plot(curWF(:, iUnit), 'k-', 'LineWidth', 3)
        ylim([-3000 3000])
        xlim([1 25])
        axis off
        
        % save
        saveFileName = sprintf('fig2c_%d%s', idChan, unitLabel{iUnit});
        print(fig2c, fullfile(dirFig, saveFileName), '-depsc')
    end
end

%% Scale bar
unitTimeSec = 0.001; % 1 ms
unitAmpMilliV = 0.050; % 50 micro volt

fig2cScale = figure; % for each unit
set(fig2cScale, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 100 300])
clf
subplot('Position', [0 0 1 1]);
line([1 1; 1 1+unitTimeSec/(1/24414)], [1000 1000; 1000+(unitAmpMilliV*(2^15)) 1000], 'Color', 'k', 'LineWidth', 3)
ylim([-3000 3000])
xlim([1 25])
axis off
 % save
 saveFileName = 'fig2c_ScaleBar';
 print(fig2cScale, fullfile(dirFig, saveFileName), '-depsc')
    
    
    
%     % Option 2: Take average waveform across sessions
%     catWF=[];
%     for iS = 1:3
%         if size(data.waveform{iS},2)>1
%             catWF = cat(3, catWF, data.waveform{iS});
%         end
%     end
%     curWF = mean(catWF, 3);
    
%     figure(figRawWF)
%     subplot(1, numChan, iChan)
%     plot(data.waveform{1}); hold on
%     plot(data.waveform{2}); hold on
%     plot(data.waveform{3}); hold on
%     title(num2str(idChan))
%     ylim([-3000 3000])
% %     xlim([0 40])
%     
%     figure(figAvgWF)
%     subplot(1, numChan, iChan)
%     plot(avgWF, 'o-')
%     % plot(data.waveform{1}); hold on
%     % plot(data.waveform{2}); hold on
%     % plot(data.waveform{3}); hold on
%     title(num2str(idChan))
%     ylim([-3000 3000])
%     xlim([0 40])
% 
% % end






%     case 'rho'
%         % Get the list of sorted files and channels
%         dPlx = dir([dirData_spike, '/*-sorted.plx']); %'/*a-sorted.plx']);
%         plxFileName={};
%         [plxFileName{1:length(dPlx),1}] = deal(dPlx.name);
%         
%         data = waveform_consistency(plxFileName, dirData_spike);        
%         
%         % save the result and the figure
%         filename = fullfile(dirResult, 'waveformComp.mat'); %sprintf(fullfile(dirResult, 'waveformComp_%d.mat'), idChan);
%         save(filename, 'data')
%         
%         % save figures
% %         waveform_visualize
%         print(gcf, fullfile(dirResult, 'waveformComp'), '-depsc') %sprintf(fullfile(dirResult, 'waveformComp_%d'), idChan), '-depsc')
% end

        