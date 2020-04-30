% xcorr_plot - plot cross-correlogram and create matrix
% 
% written by Kenji Koyano, 12/14/2015
% last modified -12/14/2015


clear all;
close all;

%% Parameters
directory.xcorr      = '/procdata/koyanok/physiology/xcorr/';

% list.tank  = {'Matcha151007';'Matcha151008';'Matcha151015';'Matcha151023';'Matcha151030'};
list.tank  = {'Matcha151119'};
list.block = {'Resting2'};

xparams.plot_range     = [-15 15]; % ms
xparams.zscore_range   = [-4 6];   % ms
xparams.peak_range     = 3;        % ms
xparams.min_spike_pair   = 1200;
xparams.min_spike_single = 300;


if exist('list','var')
    if isfield(list,'block')
        block_specified = true;
    else
        block_specified = false;
    end
end

if exist('list','var')
    if ~isfield(list,'tank')
        tmp_folderlist = dir(directory.xcorr);
        isd = [tmp_folderlist(:).isdir]; % only use folders
        list.tank = {tmp_folderlist(isd).name}';
        list.tank(ismember(list.tank,{'.','..'})) = []; % remove the . and .. from the list
        clear tmp_folderlist isd
    end
end

%% Loop for tank
for tankno = 1:length(list.tank)
    tankname = list.tank{tankno};
    disp(tankname);
    if ~block_specified
        tmp_folderlist = dir([directory.xcorr tankname filesep 'Resting*']);
        isd = [tmp_folderlist(:).isdir]; % only use folders
        list.block = {tmp_folderlist(isd).name}';
        list.block(ismember(list.block,{'.','..'})) = []; % remove the . and .. from the list
        clear tmp_folderlist isd
    end
    
    %% Loop for block
    for blockno = 1:length(list.block)
        blockname = list.block{blockno};
        disp(blockname);
        disp('loading data...')
        load([directory.xcorr tankname filesep blockname filesep 'xcorr_dat.mat']);
        disp('plotting cross-correlogram...')
        hf = figure('Position', [300 200 1100 850], 'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8],'Color','w','InvertHardcopy','off');
%         hf = figure('Position', [300 0 110*6 85*6], 'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8],'Color','w','InvertHardcopy','off');        
        ha_title = axes('Position',[0.3 0.98 0.2 0.05],'Visible','off');
        text(0,0,['Tank: ' tankname ', Block: ' blockname],'Interpreter','none','FontSize',12);
        ha_plot = axes('Position',[0.03 0.03 0.92 0.92],'XTick',[],'YTick',[]);
        xlabel('target cell','Interpreter','none','FontSize',10);
        ylabel('reference (trigger) cell','Interpreter','none','FontSize',10);
        size_plot = size(xcorr_dat.zscore_hist);
        plot_range_indx = xparams.plot_range+1-params.duration_xcorr(1);
        %% Loop for plot
        for cell1 = 1:size_plot(1)
            axes('Position',[0.05 0.95-cell1*0.9/size_plot(1) 0.9/size_plot(2)-0.005 0.9/size_plot(1)-0.005],'XTick',[],'YTick',[]);
            ylabel([num2str(xcorr_dat.list_ch(cell1)) '-' num2str(xcorr_dat.list_cluster(cell1))],'FontSize',2);
            for cell2 = 1:size_plot(2)
                if cell1 ==size_plot(1)
                    axes('Position',[0.05+(cell2-1)*0.9/size_plot(2) 0.95-cell1*0.9/size_plot(1) 0.9/size_plot(2)-0.005 0.9/size_plot(1)-0.005],'XTick',[],'YTick',[]);
                    xlabel([num2str(xcorr_dat.list_ch(cell2)) '-' num2str(xcorr_dat.list_cluster(cell2))],'FontSize',2);
                end
                ha_xcorr(cell1,cell2) = axes('Position',[0.05+(cell2-1)*0.9/size_plot(2) 0.95-cell1*0.9/size_plot(1) 0.9/size_plot(2)-0.005 0.9/size_plot(1)-0.005]);
                hp(cell1,cell2) = plot(xcorr_dat.zscore_hist{cell1,cell2},'k');
                set(gca,'XTickLabel',[],'YTickLabel',[]);
                xlim(plot_range_indx);
                ylim(xparams.zscore_range);
            end
        end
        %% Coloring
        x_indx.min_spike = xcorr_dat.n_spike>xparams.min_spike_pair & xcorr_dat.n_spike_min>xparams.min_spike_single;
        x_indx.lp = squeeze(xcorr_dat.peak_positive(:,:,1) & xcorr_dat.peak_time_positive(:,:,1)>-xparams.peak_range & ~xcorr_dat.peak_positive(:,:,2));  % left positive
        set(hp(x_indx.lp&x_indx.min_spike),'Color','r');
        x_indx.rp = squeeze(xcorr_dat.peak_positive(:,:,3) & xcorr_dat.peak_time_positive(:,:,3)<xparams.peak_range & ~xcorr_dat.peak_positive(:,:,2));   % right positive
        set(hp(x_indx.rp&x_indx.min_spike),'Color','m');
        x_indx.bp = squeeze(xcorr_dat.peak_positive(:,:,1) & xcorr_dat.peak_time_positive(:,:,1)>-xparams.peak_range & xcorr_dat.peak_positive(:,:,3) & xcorr_dat.peak_time_positive(:,:,3)<xparams.peak_range & ~xcorr_dat.peak_positive(:,:,2));  % both positive
        set(hp(x_indx.bp&x_indx.min_spike),'Color',[1 0 0.5]);
        x_indx.ln = squeeze(xcorr_dat.peak_negative(:,:,1) & xcorr_dat.peak_time_negative(:,:,1)>-xparams.peak_range & ~xcorr_dat.peak_negative(:,:,2));  % left negative
        set(hp(x_indx.ln&x_indx.min_spike),'Color','b');
        x_indx.rn = squeeze(xcorr_dat.peak_negative(:,:,3) & xcorr_dat.peak_time_negative(:,:,3)<xparams.peak_range & ~xcorr_dat.peak_negative(:,:,2));   % right negative
        set(hp(x_indx.rn&x_indx.min_spike),'Color','c');
        x_indx.bn = squeeze(xcorr_dat.peak_negative(:,:,1) & xcorr_dat.peak_time_negative(:,:,1)>-xparams.peak_range & xcorr_dat.peak_negative(:,:,3) & xcorr_dat.peak_time_negative(:,:,3)<xparams.peak_range & ~xcorr_dat.peak_negative(:,:,2));  % both negative
        set(hp(x_indx.bn&x_indx.min_spike),'Color',[0 0.5 1]);
        x_indx.cp = squeeze(xcorr_dat.peak_positive(:,:,2));  % center positive peak
        set(hp(x_indx.cp&x_indx.min_spike),'Color','g');
        x_indx.cn = squeeze(xcorr_dat.peak_negative(:,:,2));  % center negative peak
        set(hp(x_indx.cn&x_indx.min_spike),'Color',[0 1 0.5]);
        set(hp(diag(ones(length(xcorr_dat.zscore_hist),1))>0),'Color','k');   % diagonal
        x_indx.ai = abs(xcorr_dat.ai)>0.4;
        x_indx.cs = abs(xcorr_dat.cs)>20;
        set(ha_xcorr(x_indx.ai&x_indx.min_spike),'Color',[1 .95 .95]);          % asymmetry
        set(ha_xcorr(x_indx.cs&x_indx.min_spike),'Color',[.95 .95 1]);          % coupling strength
        set(ha_xcorr(x_indx.ai&x_indx.cs&x_indx.min_spike),'Color',[1 .8 .8]);  % coupling strength & asymmetry
        set(ha_xcorr(diag(ones(length(xcorr_dat.zscore_hist),1))>0),'Color','w');   % diagonal
        
        %% Save figure and data
        disp('saving figure and data...')
        print(gcf,[directory.xcorr tankname filesep blockname filesep 'xcorrelogram.png'],'-dpng');
        save([directory.xcorr tankname filesep blockname filesep 'xcorr_index.mat'],'x_indx','xparams');
    end
end
