% function [] = S_seeCorrMap()
% GUI to explore movie features related to cluster activity
% 
% by SHP 2016
% 
% 
% global STDPATH DSP DATA GH tempS curRGR curFMRITC
% global S dataBOLD CorrMapHandles 
% 
% if ~exist('S') | isempty(S)
%     disp(['Error: global S is not defined. Run the script "S_neuralRegressor"']);
%     return;
% end
% 
% if ~exist('dataBOLD') | isempty(dataBOLD)
%     disp(['Error: global dataBOLD is not defined. Run the script "S_neuralRegressor"']);
%     return;
% end

% 
% if  ~exist('CorrMapHandles') | ~isfield(CorrMapHandles,'figure') | ...
%         ~ishandle(CorrMapHandles.figure) | (isfield(CorrMapHandles, 'reset_figure') ...
%         & CorrMapHandles.reset_figure)
    
    % Initial setup for figure and panel
    ClusterMovieHandles.figure = figure(99);
    set(gcf,'Color', 'w', 'Position',[500 0 1000 650]);
    clf 
    
    % List of movies and scenes
    uicontrol('Style','Text',...
              'Units','Normalized', 'Position',[0.02 0.9 0.15 0.08],...
              'String','Movie','Parent',ClusterMovieHandles.figure,...
              'FontWeight','bold','FontSize', ...
              12);

    uicontrol('Style','Text',...
              'Position',[0.02 0.9 0.15 0.08],...
              'String','Scene','Parent',ClusterMovieHandles.figure,...
              'Units','Normalized','FontWeight','bold',...
              'ForegroundColor',[0.0 0.5 0.0],'FontSize',12);
    
    ClusterMovieHandles.SceneSelection = ...
        uicontrol('Style','Popupmenu',...
                  'Position',[0.05 0.05 0.2 0.16],...
                  'String',cells,'Parent',ClusterMovieHandles.panel,...
                  'Units','Normalized',...
                  'FontSize',12,'Callback','getCellID'); 
    
%     ClusterMovieHandles.panel =...
%         uipanel('Position',[0.05 0.66 0.3 0.3]);
   
    % List of cells and movies
    setMovie = [1 2 3];
        
    

    uicontrol('Style','Text',...
              'Position',[0.02 0.45 0.15 0.08],...
              'String','Neural','Parent',ClusterMovieHandles.panel,...
              'Units','Normalized','FontWeight','bold',...
              'ForegroundColor',[0.5 0.0 0.0],'FontSize',12);
          
    uicontrol('Style','Text',...
              'Position',[0.02 0.3 0.15 0.08],...
              'String','Eigen','Parent',ClusterMovieHandles.panel,...
              'Units','Normalized','FontWeight','bold',...
              'ForegroundColor',[0.0 0.0 0.5],'FontSize',12);
    

    ClusterMovieHandles.fMRICheck=[]; ClusterMovieHandles.NeuralCheck=[]; ClusterMovieHandles.EigenCheck=[];
    for iM=1:length(movs)
        uicontrol('Style','text','String',sprintf('M%d',movs(iM)),...
            'Units','Normalized',...
            'Position',[-0.03+0.1+ 0.8*iM/length(movs) 0.8 0.09 0.08],'Parent',...
            ClusterMovieHandles.panel);
        ClusterMovieHandles.fMRICheck(iM) = ...
            uicontrol('Style','checkbox','Value',0,'Units','Normalized',...
            'Position',[0.1+0.8*iM/length(movs) 0.6 0.03 0.06],'Parent',...
            ClusterMovieHandles.panel);
        
        ClusterMovieHandles.NeuralCheck(iM) = ...
            uicontrol('Style','checkbox','Value',0,'Units','Normalized',...
            'Position',[0.1+0.8*iM/length(movs) 0.45 0.03 0.06],'Parent',...
            ClusterMovieHandles.panel,'ForegroundColor',[0.8 1 0.8]);
        
        ClusterMovieHandles.EigenCheck(iM) = ...
            uicontrol('Style','checkbox','Value',0,'Units','Normalized',...
            'Position',[0.1+0.8*iM/length(movs) 0.3 0.03 0.06],'Parent',...
            ClusterMovieHandles.panel,'ForegroundColor',[0.8 1 0.8]);
    end
    
    

    % Compute correlation
    uicontrol('Style','Pushbutton',...
              'Position',[0.7 0.05 0.25 0.16],...
              'String','Compute','Parent',ClusterMovieHandles.panel,...
              'Units','Normalized','BackgroundColor',[0.9 0.8 0.8],...
              'FontSize',12,'Callback','computeNewMap');

    % Reset the figure
    uicontrol('Style','Pushbutton',...
              'Position',[0.4 0.05 0.2 0.13],...
              'String','reset','Parent',ClusterMovieHandles.panel,...
              'Units','Normalized','BackgroundColor',[0.8 0.8 0.9],...
              'FontSize',12,'Callback',...
              'CorrMapHandles.reset_figure=1;S_seeCorrMap');   
    ClusterMovieHandles.reset_figure = 0;
    
              
              
    % Options for computing correlation
    ClusterMovieHandles.panel2 =...
        uipanel('Position',[0.55 0.66 0.2 0.3]);
    
    % Option 1: subtract averaged fMRI signal across 9 movies    
    ClusterMovieHandles.flagSubtAvg = ...
        uicontrol('Style','checkbox','Value',0,'Units','Normalized',...
                  'Position',[0.8 0.8 0.1 0.1],'Parent',...
                  ClusterMovieHandles.panel2,'ForegroundColor',[0.8 1 0.8]);
    uicontrol('Style','Text',...
        'Position',[0.03 0.75 0.5 0.2],...
        'String','Subtract average','Parent',ClusterMovieHandles.panel2,...
        'Units','Normalized','FontWeight','bold','FontSize', 10);
              
    % Option 2: filtering 
    uicontrol('Style','Text',...
        'Position',[0.03 0.5 0.5 0.2],...
        'String','Temporal filter','Parent',ClusterMovieHandles.panel2,...
        'Units','Normalized','FontWeight','bold','FontSize', 10);
    ClusterMovieHandles.flagFiltering = ...
        uicontrol('Style','checkbox','Value',0,'Units','Normalized',...
                  'Position',[0.8 0.6 0.1 0.1],'Parent',...
                  ClusterMovieHandles.panel2,'ForegroundColor',[0.8 1 0.8]);
%     CorrMapHandles.filterType = uicontrol('Style', 'popup', 'Units', 'normalized',...
%         'Position', [0.8 0.5 0.2 0.2],...
%         'String', {'N/A', 'HPF', 'LPF'}, 'Parent', CorrMapHandles.panel2,...
%         'FontWeight','bold','FontSize', 8);
%     CorrMapHandles.filterCF = uicontrol('Style', 'popup', 'Units', 'normalized',...
%         'Position', [0.6 0.5 0.2 0.2],...
%         'String', 'N/A|0.05|0.02|0.01', 'Parent', CorrMapHandles.panel2,...
%         'FontWeight','bold','FontSize', 8);
    
%     % Option 3: make fMRI to percent change
%     CorrMapHandles.flagPercentSig = ...
%         uicontrol('Style','checkbox','Value',0,'Units','Normalized',...
%                   'Position',[0.8 0.4 0.1 0.1],'Parent',...
%                   CorrMapHandles.panel2,'ForegroundColor',[0.8 1 0.8]);
%     uicontrol('Style','Text',...
%         'Position',[0.03 0.75 0.5 0.2],...
%         'String','Subtract average','Parent',CorrMapHandles.panel2,...
%         'Units','Normalized','FontWeight','bold','FontSize', 10);


    % Option 
    uicontrol('Style','Text',...
        'Position',[0.03 0.2 0.5 0.2],...
        'String','BLP','Parent',ClusterMovieHandles.panel2,...
        'Units','Normalized','FontWeight','bold','FontSize', 10);
        ClusterMovieHandles.BLPtype = uicontrol('Style', 'popup', 'Units', 'normalized',...
        'Position', [0.6 0.05 0.3 0.3],...
        'String', 'N/A|4-12|12-25|40-60|60-150', 'Parent', ClusterMovieHandles.panel2,...
        'FontWeight','bold','FontSize', 8);
    
    
    
    % Plot the current voxel's time series and neural regressor used
    subplot('Position',[0.05 0.3 0.92 0.3]);
    dum = ones(1,125);
    [ClusterMovieHandles.both_axes, ClusterMovieHandles.fmriTC1, ClusterMovieHandles.neuralRGR1]=plotyy(dum,dum,dum,dum);
    axis(ClusterMovieHandles.both_axes, 'tight')
    legend( 'Voxel time course','Neural regressor')
    xlabel('Time (s)')
    
    
    % Plot the current voxel's filtered time series and neural regressor used
    subplot('Position',[0.05 0.05 0.92 0.2]);
    dum = ones(1,125);
    [ClusterMovieHandles.both_axes2, ClusterMovieHandles.fmriTC2, ClusterMovieHandles.neuralRGR2]=plotyy(dum,dum,dum,dum);
    axis(ClusterMovieHandles.both_axes2, 'tight')
    legend( 'Voxel time course','Neural regressor')
    xlabel('Time (s)')
    
    ClusterMovieHandles.xcorrmap  = subplot('Position',[0.8 0.66 0.15 0.3]);
    ClusterMovieHandles.xcorrdata = plot([0 0],[0 1]);
    hold on
    plot([-30 30],[0 0],'k')
    set(gca,'YLim',[-0.6 0.6]);
    plot([0 0],get(gca,'YLim'),'k');
    set(gca,'XLim',[-30 30])
    
    ClusterMovieHandles.corrval = uicontrol('Style','text','String','r=test',...
              'Position',[0.85 0.62 0.12 0.03],'Units','Normalized','FontWeight',...
                  'Bold','FontSize',12);
              
    ClusterMovieHandles.visROI = uicontrol('Style','text','String','VisROI',...
              'Position',[0.05 0.62 0.12 0.03],'Units','Normalized',...
              'ForegroundColor',[0.0 0.5 0.0],'FontWeight','Bold','FontSize',12);
              
    ClusterMovieHandles.faceROI = uicontrol('Style','text','String','FaceROI',...
              'Position',[0.25 0.62 0.12 0.03],'Units','Normalized',...
              'ForegroundColor',[0.5 0.0 0.0],'FontWeight','Bold','FontSize',12);
    return;
end

updateVolumeTC;

% Plot the current voxel's time series and neural regressor used
z_timecourse = (DSP.timecourse-nanmean(DSP.timecourse))./nanstd(DSP.timecourse);
z_currgr     = (DATA.stim.curstimmodel-nanmean(DATA.stim.curstimmodel))./nanstd(DATA.stim.curstimmodel);

set(CorrMapHandles.fmriTC1,'XData',DATA.ftimes); %set(get(CorrMapHandles.both_axes(1),'Children'),'XData',DATA.ftimes);
set(CorrMapHandles.fmriTC1,'YData',z_timecourse); %set(get(CorrMapHandles.both_axes(1),'Children'),'YData',z_timecourse);
set(CorrMapHandles.both_axes(1),'YLim',[-4 4]);
set(CorrMapHandles.both_axes(1),'XLim',[DATA.ftimes(1) DATA.ftimes(end)]);

set(CorrMapHandles.neuralRGR1,'XData',DATA.ftimes); %set(get(CorrMapHandles.both_axes(2),'Children'),'XData',DATA.ftimes);
set(CorrMapHandles.neuralRGR1,'YData',z_currgr); %set(get(CorrMapHandles.both_axes(2),'Children'),'YData',z_currgr);
set(CorrMapHandles.both_axes(2),'YLim',[-4 4]);
set(CorrMapHandles.both_axes(2),'XLim',[DATA.ftimes(1) DATA.ftimes(end)]);


% Find the actual correlation value
arr = DSP.proc.scalarmap_3d;
voxi = DSP.or.activevoxel;
set(CorrMapHandles.corrval,'String',sprintf('r = %.2f',arr(voxi(1),voxi(2),voxi(3))));

% Display ROI if the active voxel belongs to any
[flagVisROI, flagFaceROI, nameROI] = determineROIs(voxi);
if flagVisROI
    set(CorrMapHandles.visROI,'String', nameROI); %sprintf('r = %.2f',arr(voxi(1),voxi(2),voxi(3))));
    set(CorrMapHandles.faceROI,'String','FaceROI'); %
end
if flagFaceROI
    set(CorrMapHandles.visROI,'String', 'VisROI'); %sprintf('r = %.2f',arr(voxi(1),voxi(2),voxi(3))));
    set(CorrMapHandles.faceROI,'String', nameROI); %
end
if ~flagVisROI & ~flagFaceROI
    set(CorrMapHandles.visROI,'String', 'VisROI'); %
    set(CorrMapHandles.faceROI,'String', 'FaceROI'); %
end

% Plot the filtered current voxel's time series and neural regressor used
z_timecourse_filtered = (curFMRITC(voxi(1), voxi(2), voxi(3), :)...
    -nanmean(curFMRITC(voxi(1), voxi(2), voxi(3), :)))./nanstd(curFMRITC(voxi(1), voxi(2), voxi(3), :));
% mion
if DSP.mion, z_timecourse_filtered = -z_timecourse_filtered;end

z_currgr_filtered     = (curRGR-nanmean(curRGR))./nanstd(curRGR);

set(CorrMapHandles.fmriTC2,'XData',DATA.ftimes); %set(get(CorrMapHandles.both_axes2(1),'Children'),'XData',DATA.ftimes);
set(CorrMapHandles.fmriTC2,'YData',z_timecourse_filtered); %set(get(CorrMapHandles.both_axes2(1),'Children'),'YData',z_timecourse_filtered);
set(CorrMapHandles.both_axes2(1),'YLim',[-4 4]);
set(CorrMapHandles.both_axes2(1),'XLim',[DATA.ftimes(1) DATA.ftimes(end)]);

set(CorrMapHandles.neuralRGR2,'XData',DATA.ftimes); %set(get(CorrMapHandles.both_axes2(2),'Children'),'XData',DATA.ftimes);
set(CorrMapHandles.neuralRGR2,'YData',z_currgr_filtered); %set(get(CorrMapHandles.both_axes2(2),'Children'),'YData',z_currgr_filtered);
set(CorrMapHandles.both_axes2(2),'YLim',[-4 4]);
set(CorrMapHandles.both_axes2(2),'XLim',[DATA.ftimes(1) DATA.ftimes(end)]);

% for now just drawing all this...maybe convert to handles later

% Cross-correlogram
timecourse = squeeze(curFMRITC(voxi(1),voxi(2),voxi(3),:));
mntc = nanmean(timecourse);
timecourse(isnan(timecourse)) = mntc;

% -1 below for mion
[xc,lags] = xcorr(zscore(curRGR),-1*zscore(timecourse),20,'coeff');
[TR,NR] = getTRandNR();
lags = lags*TR;
set(CorrMapHandles.xcorrdata,'XData',lags,'YData',xc)


%axis(CorrMapHandles.both_axes, 'tight')






