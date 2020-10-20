% genFig_waveform.m
%
% Plot the waveforms of units from the same channels
% For a given set of units, 1) Read the .plx files for each unit and 
% 2) Collect the waveform data then 3) Plot
% 2015/12/16 SHP

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
dirNeuralData = fullfile(dirProcdata, '/parksh/_macaque', nameSubject);
dirData_waveform = fullfile(dirNeuralData, '/_sortedSpikes/_waveformCompare');
dirData_spike = fullfile(dirNeuralData, '_sortedSpikes');
 

% toolkit for reading .plx file
addpath(fullfile(dirProjects, '/parksh/_toolbox/Waveform_Consistency-V1.1/readPLXFileC'))

%% Which units?
listChan = [82 85 86 89 91 96 122]'; 
numChan = length(listChan);

%% Read files and draw the waveform (for all channels)
% setSessions = {'toroid20120617b', 'toroid20120618b', 'toroid20120619a'}; % for movie 123
figAvgWF = figure;
set(figAvgWF, 'Color', 'w', 'PaperPositionMode', 'auto')
figRawWF = figure;
set(figRawWF, 'Color', 'w', 'PaperPositionMode', 'auto')
for iChan = 1:numChan
    plxFileName = {};
    idChan = listChan(iChan);
    filename_wf = fullfile(dirData_waveform, sprintf('waveformComp_%d', idChan));
    load(filename_wf)
    % data.merged_channel_units
    
    catWF=[];
    for iS = 1:3
        if size(data.waveform{iS},2)>1
            catWF = cat(3, catWF, data.waveform{iS});
        end
    end
    avgWF = mean(catWF, 3);
    
    figure(figRawWF)
    subplot(1, numChan, iChan)
    plot(data.waveform{1}); hold on
    plot(data.waveform{2}); hold on
    plot(data.waveform{3}); hold on
    title(num2str(idChan))
    ylim([-3000 3000])
%     xlim([0 40])
    
    figure(figAvgWF)
    subplot(1, numChan, iChan)
    plot(avgWF, 'o-')
    % plot(data.waveform{1}); hold on
    % plot(data.waveform{2}); hold on
    % plot(data.waveform{3}); hold on
    title(num2str(idChan))
    ylim([-3000 3000])
    xlim([0 40])

end






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

        