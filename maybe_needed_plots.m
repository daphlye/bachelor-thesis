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