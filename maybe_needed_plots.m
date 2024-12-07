%% particle distribution time series
% time series ele

nPanels = 8;
h = irf_plot(nPanels);
   % tmp.nRows = 4; tmp.nCols = 2;
   % for ii = 1:tmp.nRows*tmp.nCols; h(ii) = subplot(tmp.nRows,tmp.nCols,ii); end
   % tmp.isub = 1;
    %
%h = irf_plot(tmp.nPanels);
c_eval('event.ePDistN = event.ePDist?; event.raw.BN_gsm = event.raw.B?_gsm;',Param.ic)
% B
if 1
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('B');
  irf_plot(hca,event.raw.BN_gsm)
  hold(hca,'on')
  irf_plot(hca,event.raw.BN_gsm.abs)
  hca.YLabel.String = 'B (nT)';
  title(hca,irf_ssub('MMS? electrons', Param.ic))
end
% J*r
%hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
hca = irf_panel('e omni');
irf_plot(hca, event.jr_mag.all, 'o', 'MarkerSize', 3);
ylabel(hca, '$\vec{\mathbf{J}} \cdot \hat{\mathbf{r}}$', 'Interpreter', 'latex');

% e omni
if 1
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('e omni 64 energy channels'); %#ok<UNRCH>
  [hca, hcb]=irf_spectrogram(hca,event.ePDistN.omni.specrec);
  colormap('jet');
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
end
% e pitchangles lower time resolution
if 1
  elim = [20 200];
  hca = irf_panel('e pitchangles 64 energy channels');
  irf_spectrogram(hca,event.ePDistN.e64.elim(elim).pitchangles(event.raw.BN_gsm,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
% e  pa mid
if 1
  elim = [200 2000];
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('e pitchangles mid');
  irf_spectrogram(hca,event.ePDistN.e64.elim(elim).pitchangles(event.raw.BN_gsm,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
% e pa high
if 0
  elim = [400 20000]; %#ok<UNRCH>
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('e pitchangles high');
  irf_spectrogram(hca,event.ePDistN.elim(elim).pitchangles(event.raw.BN_gsm,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
% e spectogram parallel
if 1
  palim = [0 15];
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
    hca = irf_panel('e spectrogram parallel');
  irf_spectrogram(hca,event.ePDistN.pitchangles(event.raw.BN_gsm,[0 15]).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
% e spectogram perpendicular
if 1
  palim = [75 105];
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
    hca = irf_panel('e spectrogram perpendicular');
  irf_spectrogram(hca,event.ePDistN.pitchangles(event.raw.BN_gsm,palim).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
% e spectogram anti parallel
if 1
  palim = [165 180];
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('e spectrogram anti-parallel');
  irf_spectrogram(hca,event.ePDistN.pitchangles(event.raw.BN_gsm,palim).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end

irf_zoom(h,'x',tint_beob)
irf_plot_axis_align
tmp.filename = fullfile('events/', event.tint_string,'electron_distributions');
irf_print_fig(tmp.filename,'png');
close(gcf);
%% time series ion
nPanels = 8;
h = irf_plot(nPanels);
   % tmp.nRows = 4; tmp.nCols = 2;
   % for ii = 1:tmp.nRows*tmp.nCols; h(ii) = subplot(tmp.nRows,tmp.nCols,ii); end
   % tmp.isub = 1;
    %
%h = irf_plot(tmp.nPanels);
c_eval('event.iPDistN = event.iPDist?;',Param.ic)
% B
if 1
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('B');
  irf_plot(hca,event.raw.BN_gsm)
  hold(hca,'on')
  irf_plot(hca,event.raw.BN_gsm.abs)
  hca.YLabel.String = 'B (nT)';
  title(hca,irf_ssub('MMS? Ions', Param.ic))
end
% J*r
%hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
hca = irf_panel('p omni');
irf_plot(hca, event.jr_mag.all, 'o', 'MarkerSize', 3);
ylabel(hca, '$\vec{\mathbf{J}} \cdot \hat{\mathbf{r}}$', 'Interpreter', 'latex');

% p omni
if 1
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('p omni 64 energy channels'); %#ok<UNRCH>
  irf_spectrogram(hca,event.iPDistN.omni.specrec)
  colormap('jet');
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
end

% e pitchangles lower time resolution
if 1
  elim = [20 200];
  hca = irf_panel('p pitchangles 64 energy channels');
  irf_spectrogram(hca,event.iPDistN.e64.elim(elim).pitchangles(event.raw.BN_gsm,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_p<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
% e  pa mid
if 1
  elim = [200 2000];
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('p pitchangles mid');
  irf_spectrogram(hca,event.iPDistN.e64.elim(elim).pitchangles(event.raw.BN_gsm,18).specrec('pa'));
  hca.YTick = [0 45 90 135];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_p<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end

% e spectogram parallel
if 1
  palim = [0 15];
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
    hca = irf_panel('p spectrogram parallel');
  irf_spectrogram(hca,event.iPDistN.pitchangles(event.raw.BN_gsm,[0 15]).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
% e spectogram perpendicular
if 1
  palim = [75 105];
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
    hca = irf_panel('p spectrogram perpendicular');
  irf_spectrogram(hca,event.iPDistN.pitchangles(event.raw.BN_gsm,palim).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
% e spectogram anti parallel
if 1
  palim = [165 180];
  %hca = h(tmp.isub); tmp.isub = tmp.isub + 1;
  hca = irf_panel('p spectrogram anti-parallel');
  irf_spectrogram(hca,event.iPDistN.pitchangles(event.raw.BN_gsm,palim).specrec,'log')
  hca.YScale = 'log';
  hca.YTick = 10.^[1 2 3 4];
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end

irf_zoom(h,'x',tint_beob)
irf_plot_axis_align;
tmp.filename = fullfile('events/', event.tint_string,'ion_distributions');
irf_print_fig(tmp.filename,'png');
close(gcf);

%% velocity distributions electrons
if 0
    %stuff to plot
    event.vd.t0 = tint_beob(1); %2015-12-02T01:14:54.320Z
    % tInd is timeindex, where t0 in downloaded tints in
    % ePDist? begins:
    c_eval('event.vd.tInd = find(abs(event.ePDist?.time-event.vd.t0)==min(abs(event.ePDist?.time-event.vd.t0)));', Param.ic); 
    c_eval('event.vd.B0 = event.raw.B?_gsm.resample(event.ePDist?.time).data;',Param.ic);
    c_eval('event.vd.E0 = event.dslE?.resample(event.ePDist?.time).data;',Param.ic);
    c_eval('event.scpot = event.scPot?.resample(event.ePDist?.time);',Param.ic);
    c_eval('event.ePitchN = event.ePDist?.pitchangles(event.raw.B?_gsm,[17]);',Param.ic);
    %c_eval('event.iPDist?_res = event.iPDist?.resample(event.ePDist?.time).data;',Param.ic);
    % Is already at same sampling rate since measured by same instrument???
    c_eval('event.iPDistN = event.iPDist?;',Param.ic);
    c_eval('event.ePDistN = event.ePDist?;',Param.ic);
    event.hatB0 = double(irf_norm(event.vd.B0));
    event.hatE0tmp = double(irf_norm(event.vd.E0));
end
%% input parameters for projection plot
Param.vlim = 12*1e3; % x and ylim
Param.vlim_i = 70*1e1; % x and ylim
Param.elevlim = 15; % angle over plane to include in slice
Param.strCMap = 'jet'; % colormap
Param.projclim = [0 5]; % colorbar limit
Param.projclim_i = [0 10]; % still colorbar limit but for the ions
%% E ExB B direction velocity distribution
if 0
    tmp.filename_dist = fullfile('events/', event.tint_string, 'desdist_e_EBExB');
    if ~exist(tmp.filename_dist, 'dir')
        mkdir(tmp.filename_dist);
    end
    
    Param.step = 1;  % Tipps: step = 30 => 1/s
    for i = Param.selFr 
      
      tmp.idx = event.vd.tInd + (i-1)*Param.step; % index
      tmp.time = event.ePDistN.time(tmp.idx);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % v dist in E ExB B direction
      if 1
          event.hatExB0 = cross(event.hatE0tmp(tmp.idx,:),event.hatB0(tmp.idx,:));
          event.hatE0 = cross(event.hatB0(tmp.idx,:),event.hatExB0);
          %define Axis, since its in the direction of E etc its getting calculated
          %new for every time! 
          x = event.hatE0;
          y = event.hatExB0;
          z = event.hatB0(tmp.idx,:);
          
          % Initialize figure
          figHandle = figure('Visible', 'off');
          % Alt 1: 2x3 plots
          %nRows = 2; nCols = 3;
          %for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end
          % Alt 2: 3 + 2 plots
          for ii = 1:3; h(ii) = subplot(2,3,ii); end
          for ii = 4:5; h(ii) = subplot(2,2,ii-1); end
          isub = 1;
          
          hca = h(isub); isub = isub + 1;
          xyz = [x; y; z * (-1)];
          vlabels = {'v_E','v_{ExB}','v_B'};
          mms.plot_projection(hca,event.ePDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim,'clim',Param.projclim,'scpot',event.scpot,'vlabel',vlabels);
          
          hca = h(isub); isub = isub + 1;
          xyz = [y; z; x * (-1)];
          vlabels = {'v_{ExB}','v_B','v_E'};
          mms.plot_projection(hca,event.ePDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim,'clim',Param.projclim,'scpot',event.scpot,'vlabel',vlabels);
          
          hca = h(isub); isub = isub + 1;
          xyz = [z; x; y * (-1)];
          vlabels = {'v_B','v_E','v_{ExB}'};
          mms.plot_projection(hca,event.ePDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim,'clim',Param.projclim,'scpot',event.scpot,'vlabel',vlabels);
          
          if 0
            hca = h(isub); isub = isub + 1; %#ok<UNRCH>
            mms.plot_skymap(hca,event.ePDistN,'tint',tmp.time,'energy',150,'flat');
            
            hca = h(isub); isub = isub + 1;
            mms.plot_skymap(hca,event.ePDistN,'tint',tmp.time,'energy',150,'flat','log');
            %hca.CLim = Param.projclim;
            
            hca = h(isub); isub = isub + 1;
            mms.plot_skymap(hca,event.ePDistN,'tint',tmp.time,'energy',150,'vectors',{hatB0,'B'},'log');
          end
          
          hca = h(isub); isub = isub + 1;
          plot(hca,event.ePitchN.depend{1}(tmp.idx,:),event.ePitchN.data(tmp.idx,:,1),...
            event.ePitchN.depend{1}(tmp.idx,:),event.ePitchN.data(tmp.idx,:,9),...
            event.ePitchN.depend{1}(tmp.idx,:),event.ePitchN.data(tmp.idx,:,17));
          hca.XScale = 'log';
          hca.YScale = 'log';
          hca.YLabel.String = ['f_e (' event.ePitchN.units ')'];
          hca.XLabel.String = 'E_e (eV)';
          hca.XLim = [1e1 1e3];
          hca.YLim = [1e-31 2e-26];
          hleg = irf_legend(hca,{'0';'90';'180'},[0.98 0.98]);
          title(hca,irf_ssub('MMS?',Param.ic))
          
          hca = h(isub); isub = isub + 1;
          plot(hca,event.ePitchN.depend{2},squeeze(event.ePitchN.data(tmp.idx,:,:)));
          hca.YScale = 'log';
          hca.YLabel.String = ['f_e (' event.ePitchN.units ')'];
          hca.XLabel.String = '\theta (deg)';
          hca.XLim = [0 180];
          hca.YLim = [1e-31 2e-26];
          hca.XTick = [0 45 90 135 180];
          
          tmp.filename = fullfile('events/', event.tint_string, 'desdist_e_EBExB/', [irf_fname(tmp.time,4) '.png']);
          saveas(figHandle, tmp.filename);
          close(figHandle);
      end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v dist in x y z direction 
%% electron and ion velocity distributions
% input parameters for projection plot
Param.vlim = 8*1e3; % x and ylim
Param.vlim_i = 60*1e1; % x and ylim
Param.elevlim = 15; % angle over plane to include in slice
Param.strCMap = 'jet'; % colormap
Param.projclim = [0 6]; % colorbar limit
Param.projclim_i = [0 10]; % still colorbar limit but for the ions

% create folder to safe file in
tmp.filename_dist = fullfile('events/', event.tint_string, 'desdist_eandi_xyz');
if ~exist(tmp.filename_dist, 'dir')
    mkdir(tmp.filename_dist);
end

Param.selFr = randperm(100, 3);  % selected frames that will be plotted, random for each event, but same for all plots in that event
for i = Param.selFr 
    tmp.idx = event.vd.tInd + i; %tInd is first measurement point of tint
    tmp.time = event.ePDistN.time(tmp.idx);
    if 1
          x = [1,0,0];
          y = [0,1,0];
          z = [0,0,1];
          
          % Initialize figure
          figHandle = figure('Visible', 'off');
          for ii = 1:3; h(ii) = subplot(2,3,ii); end
          for ii = 4:6; h(ii) = subplot(2,3,ii); end
          isub = 1;
          %title('electrons');
          hca = h(isub); isub = isub + 1;
          colormap('jet');
          xyz = [x; y; z ];% * (-1)
          vlabels = {'v_x','v_y','v_z'};
          mms.plot_projection(hca,event.ePDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim,'clim',Param.projclim,'scpot',event.scpot,'vlabel',vlabels);
          
          hca = h(isub); isub = isub + 1;
          xyz = [y; z; x ]; % x * (-1)
          vlabels = {'v_y','v_z','v_x'};
          mms.plot_projection(hca,event.ePDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim,'clim',Param.projclim,'scpot',event.scpot,'vlabel',vlabels);
          
          hca = h(isub); isub = isub + 1;
          xyz = [z; x; y ];%* (-1)
          vlabels = {'v_z','v_x','v_y'};
          mms.plot_projection(hca,event.ePDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim,'clim',Param.projclim,'scpot',event.scpot,'vlabel',vlabels);
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% now ions
          hca = h(isub); isub = isub + 1;
          colormap('jet');
          xyz = [x; y; z ];% * (-1)
          vlabels = {'v_x','v_y','v_z'};
          mms.plot_projection(hca,event.iPDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim_i,'clim',Param.projclim_i,'scpot',event.scpot,'vlabel',vlabels);
          hca = h(isub); isub = isub + 1;
          xyz = [y; z; x ]; % x * (-1)
          vlabels = {'v_y','v_z','v_x'};
          mms.plot_projection(hca,event.iPDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim_i,'clim',Param.projclim_i,'scpot',event.scpot,'vlabel',vlabels);
          
          hca = h(isub); isub = isub + 1;
          xyz = [z; x; y ];%* (-1)
          vlabels = {'v_z','v_x','v_y'};
          mms.plot_projection(hca,event.iPDistN.convertto('s^3/km^6'),'tint',tmp.time,'xyz',xyz,'elevationlim',Param.elevlim,'vlim',Param.vlim_i,'clim',Param.projclim_i,'scpot',event.scpot,'vlabel',vlabels);

          tmp.filename = fullfile('events/', event.tint_string, 'desdist_eandi_xyz/', [irf_fname(tmp.time,4) '.png']);
          saveas(figHandle, tmp.filename);
          close(figHandle);
    end
end