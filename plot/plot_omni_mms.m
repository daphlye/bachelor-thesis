%% plot_omni_mms;
%% plot OMNI
h = irf_plot(4,'newfigure');

hca = irf_panel('M_{MS, OMNI}');
irf_plot(hca,omni.Ms);
ylabel(hca,{'M_{MS}','(#)'},'Interpreter','tex');
irf_zoom(hca,'x',tint);
irf_zoom(hca,'y',[0 5]);
irf_legend(hca,'(a)',[0.99 0.98],'color','k')

hca = irf_panel('B_{OMNI}');
irf_plot(hca,omni.Bxyz);
ylabel(hca,{'B_{GSM}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
%irf_zoom(hca,'x',tint);
irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca = irf_panel('v_{OMNI}');
irf_plot(hca,omni.Vxyz);
ylabel(hca,{'v_{GSM}','(km/s)'},'Interpreter','tex');
irf_legend(hca,{'v_{x}','v_{y}','v_{z}'},[0.88 0.10])
%irf_zoom(hca,'x',tint);
irf_legend(hca,'(c)',[0.99 0.98],'color','k')

hca = irf_panel('ni_{OMNI}'); % set(hca,'ColorOrder',mmsColors)
irf_plot(hca,omni.np);
ylabel(hca,{'n_{p}','(#/cm^{3})'},'Interpreter','tex'); % Units confirmed (SPDF)
%set(hca,'yscale','log');
irf_zoom(hca,'y',[0 10]);
irf_legend(hca,'(d)',[0.99 0.98],'color','k')
irf_plot_axis_align(1,h(1:4))
irf_zoom(h(1:4),'x',tint);
%saveas(h, 'plots/tint1/omni');  % Save the current figure as a .fig file
tmp.filename = fullfile('events/plots/','omni/');
if ~exist(tmp.filename, 'dir')
    mkdir(tmp.filename);
end
tmp.filename = fullfile('events/plots/omni/',event.tint_string);
irf_print_fig(tmp.filename,'png');
close(gcf);

%% plot MMS
if 0
    h = irf_plot(4,'newfigure');
    
    hca = irf_panel('B_{MMS1}');
    irf_plot(hca,mmsd.B1_gsm);
    ylabel(hca,{'B_{GSM}','(nT)'},'Interpreter','tex');
    irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
    %irf_zoom(hca,'y',[-70 70]);
    irf_legend(hca,'(a)',[0.99 0.98],'color','k')
    %irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.99 0.1],'color','cluster')
    
    hca = irf_panel('V_{MMS1}');
    irf_plot(hca,mmsd.Vi1_gsm);
    ylabel(hca,{'v_{GSM}','(km/s)'},'Interpreter','tex');
    irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
    %irf_zoom(hca,'y',[-70 70]);
    irf_legend(hca,'(b)',[0.99 0.98],'color','k')
    
    mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
    hca = irf_panel('n_{p, MMS}'); set(hca,'ColorOrder',mmsColors)
    irf_pl_tx(hca,'mmsd.Ni?');
    ylabel(hca,{'n_p','(#/cm^{3})'},'Interpreter','tex');
    set(hca,'yscale','log'); % ,'log'
    %irf_zoom(hca,'y',[1e-4 10]); % [1e-4 10]
    irf_legend(hca,'(c)',[0.99 0.98],'color','k')
    irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.99 0.1],'color','cluster')
    
    hca = irf_panel('J_{MMS}');
    irf_plot(hca,mmsd.j*1e9);
    ylabel(hca,{'J_{GSM, cm}','(nA/m^{2})'},'Interpreter','tex');
    irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
    irf_legend(hca,'(d)',[0.99 0.98],'color','k')
    
    if 0
        hca = irf_panel('divovercurl');
        irf_plot(hca,mmsd.divovercurl);
        ylabel(hca,{'|\nabla . B|','|\nabla \times B|'},'Interpreter','tex');
        irf_legend(hca,'(d)',[0.99 0.98],'color','k')
    end
    if 0
        hca = irf_panel('J_{MMS_pmm}');
        irf_plot(hca,mmsd.J_pmm*1e9);
        ylabel(hca,{'J_{GSM, pmm}','(nA/m^{2})'},'Interpreter','tex');
        irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
        irf_legend(hca,'(d)',[0.99 0.98],'color','k')
    end
    
    %title(h(1),'MMS - Current density and fields');
    
    irf_plot_axis_align(1,h(1:4))
    irf_zoom(h(1:4),'x',tint);
    
    tmp.filename = fullfile('events/plots/','mms/');
    if ~exist(tmp.filename, 'dir')
        mkdir(tmp.filename);
    end
    tmp.filename = fullfile('events/plots/mms/', event.tint_string);
    irf_print_fig(tmp.filename,'png');
    close(gcf);
end