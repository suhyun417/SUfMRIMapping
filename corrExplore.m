function corrExplore()

global DATA GH DSP EVT
events

[TR,NR] = getTRandNR();

% DSP.mode.correxploremode = get(GH.but.correxplore,'Value');

GH.tce.ax = [];

vol1 = DSP.proc.fmri_tc_3d;
ntr = size(vol1,4);
% DATA.tcexplore.allindxs = [DSP.TR_Discard:ntr];


%%%
  % If we are supposed to be here in the first place
  %
  if ~DSP.mode.correxploremode
      set(GH.figs.correxplore,'Visible','off');
      updateFigBar;
      return;
  else
      figure(GH.figs.correxplore);
      clf
      set(GH.figs.correxplore,'Visible','on');
      %	DSP.timecourses = [];
      
      switch DSP.lb.ori
          case 'hor'
              w = DSP.proc.params3d.dims(1);
              h = DSP.proc.params3d.dims(2);
          case 'cor'
              w = DSP.proc.params3d.dims(1);
              h = DSP.proc.params3d.dims(3);
          case 'sag'
              w = DSP.proc.params3d.dims(2);
              h = DSP.proc.params3d.dims(3);
      end
      
      
      
      if isfield(DATA,'dgzs') & ~isempty(DATA.dgzs) & isstruct(DATA.dgzs)
          
          nfiles = length(DATA.dgzs);
          yl = 8;
          
          for f=1:nfiles
              % get the eye movement data
              nobs   = length(DATA.dgzs(f).obs_times);
              [hems,vems,times] = getEVTEMs(DATA.dgzs(f));
              endobst = DATA.dgzs(f).obs_times(nobs)+length(hems{nobs});
              DATA.heyetrace{f} = nan(1,endobst);
              DATA.veyetrace{f} = nan(1,endobst);
              for i=1:nobs
                  obsstart  = DATA.dgzs(f).obs_times(i)+1;
                  obsend    = obsstart+length(hems{i})-1;
                  DATA.heyetrace{f}([obsstart:obsend]) = hems{i};
                  DATA.veyetrace{f}([obsstart:obsend]) = vems{i};
              end
              DATA.heyetrace{f}(abs(DATA.heyetrace{f})>yl) = NaN;
              DATA.veyetrace{f}(abs(DATA.veyetrace{f})>yl) = NaN;
          end
          
          % plot the eye movement data H
          GH.tce.ax.emh = subplot('Position',[0.2 0.28 0.63 0.18]);
          for f=1:nfiles
              tvals = 5*[1:length(DATA.heyetrace{f})]/1000;
              plot(tvals,DATA.heyetrace{f},'Color',getTraceCol(f),'LineWidth',0.1);
              hold on
          end
          
          axis tight
          set(gca,'YLim',[-yl yl],'Color',get(gcf,'Color'),'XTickLabel',[]);
          ylabel('H Deg');
          GH.tce.ax.emv = subplot('Position',[0.2 0.06 0.63 0.18]);
          % plot the eye movement data V
          for f=1:nfiles
              tvals = 5*[1:length(DATA.veyetrace{f})]/1000;
              plot(tvals,DATA.veyetrace{f},'Color',getTraceCol(f),'LineWidth',0.1);
              hold on
          end
          axis tight
          set(gca,'YLim',[-yl yl],'Color',get(gcf,'Color'));
          ylabel('V Deg');
      end
      
      
      %%%%%%%%%%%%%%
      %
      % plot the voxel time course data in wide form
      %
      %
      
      GH.tce.ax.data = subplot('Position',[0.2 0.52 0.63 0.4]);
      GH.tce.ax.narrow = subplot('Position',[0.03 0.06 0.13 0.87]);
      
      if isfield(DSP.lb,'tcs') & ~isempty(DSP.lb.tcs)
          nfiles = size(DSP.lb.tcs,4);
      else
          nfiles = 1;
      end
      for i=1:nfiles
          if nfiles == 1
              ftc = DSP.timecourse;
          else
              ftc = DSP.timecourses(i,:);
          end
          axes(GH.tce.ax.data);
          hold on
          GH.tce.ax.tc(i) = plot(DATA.ftimes,ftc,'Color',getTraceCol(i),'LineWidth',2);
          axes(GH.tce.ax.narrow);
          hold on
          GH.tce.ax.tc_narrow(i) = plot(DATA.ftimes,ftc,'Color',getTraceCol(i),'LineWidth',1);
          run{i} = sprintf('run %d',i);
      end
      
      if strcmp(DATA.stim.curstattype,'rgr')
          %  Draw the stimulus model
          %
          axes(GH.tce.ax.data);
          GH.tce.ax.tc(2) = plot(DATA.ftimes,zeros(1,length(DATA.ftimes)),'Color','k','LineWidth',2);
          
          axes(GH.tce.ax.narrow);
          GH.tce.ax.tc_narrow(2) = plot(DATA.ftimes,zeros(1,length(DATA.ftimes)),'Color','k','LineWidth',2);
          run{2} = sprintf('RGR');
      end
      
      axes(GH.tce.ax.data);
      GH.tce.stimdur = rectangle('Position',[0 0 0.1 0.1],'FaceColor',[0.4 0.1 0.1]);
      set(GH.tce.stimdur,'Visible','off');
      set(gca,'Color',get(gcf,'Color'),'XTickLabel',[]);
      set(gca,'YLim',[-DSP.tcscale DSP.tcscale]);
      axes(GH.tce.ax.data);
      legend(run)
      
      axes(GH.tce.ax.narrow);
      GH.tce.stimdur_narrow = rectangle('Position',[0 0 0.1 0.1],'FaceColor',[0.4 0.1 0.1]);
      set(GH.tce.stimdur_narrow,'Visible','off');
      
      set(gca,'Color',get(gcf,'Color'),'XTickLabel',[]);
      set(gca,'YLim',[-DSP.tcscale DSP.tcscale]);
      if strcmp(DATA.stim.curstattype,'rgr')
          GH.tce.ax.tc_narrow(2) = plot(DATA.ftimes,zeros(1,length(DATA.ftimes)),'Color','k','LineWidth',2);
      end
  end
  
  
  
  
  
  
  
  