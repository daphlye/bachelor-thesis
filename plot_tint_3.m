function [plt] = plot_tint_3(omni, mmsd, event, tint_beob, tint, Param)
%% mms and omni raw for whole extracted time 
% Param.ic ist das sc das wir untersuchen! also als Parameter in ide Funktion mit
% eingeben!
% create filepath
event.tint_string = irf_fname(tint_beob);
tmp.filename = fullfile('events/', event.tint_string);
if ~exist(tmp.filename, 'dir')
    mkdir(tmp.filename);
end
% set linewidth and fontsize 
set(groot, 'DefaultLineLineWidth', 2); % Change 1.5 to desired width
set(groot, 'DefaultAxesFontSize', 20);

%plt = 0; if time, we can make this a finction that tests if all plots
%got created so that code can stop if not.

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
tmp.filename = fullfile('events/', event.tint_string,'omni');
irf_print_fig(tmp.filename,'png');
close(gcf);

%% plot MMS
h = irf_plot(5,'newfigure');

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

hca = irf_panel('J_{MMS_pmm}');
irf_plot(hca,mmsd.J_pmm*1e9);
ylabel(hca,{'J_{GSM, pmm}','(nA/m^{2})'},'Interpreter','tex');
irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
irf_legend(hca,'(d)',[0.99 0.98],'color','k')

%title(h(1),'MMS - Current density and fields');

irf_plot_axis_align(1,h(1:5))
irf_zoom(h(1:5),'x',tint);
tmp.filename = fullfile('events/', event.tint_string,'mms');
irf_print_fig(tmp.filename,'png')
close(gcf);
%% sc position in R_E
figure('Visible','off');
[tmp.X,tmp.Y,tmp.Z] = sphere;
hold on;
plot3(0,0,0,'o');
hold on;
plot3(mmsd.R1_gsm.data(:,1),mmsd.R1_gsm.data(:,2),mmsd.R1_gsm.data(:,3));
plot3(mmsd.R2_gsm.data(:,1),mmsd.R2_gsm.data(:,2),mmsd.R2_gsm.data(:,3));
plot3(mmsd.R3_gsm.data(:,1),mmsd.R3_gsm.data(:,2),mmsd.R3_gsm.data(:,3));
plot3(mmsd.R4_gsm.data(:,1),mmsd.R4_gsm.data(:,2),mmsd.R4_gsm.data(:,3));

% Sphere
[tmp.p,tmp.q,tmp.r]=sphere;
tmp.p=tmp.p*6378;
tmp.q=tmp.q*6378;
tmp.r=tmp.r*6357;
surf(tmp.p,tmp.q,tmp.r);
hold off;
axis equal;
xlabel('X/km');
ylabel('Y/km');
zlabel('Z/km');
grid on;
close;
% TASK: plot Gabriellas position show code

%% all relevant data of that certain event

h = irf_figure(1,4,'reset');

hca = irf_panel('B_{MMS1}');
irf_plot(hca,event.raw.B1_gsm);
ylabel(hca,{'B_{GSM, MMS1}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])

hca = irf_panel('R_{MMS1}');
irf_plot(hca,TSeries(event.raw.R1_gsm.time, (event.raw.R1_gsm.data(:,1:3)/6378)));
ylabel(hca,{'R_{GSM, MMS1}','(R_E)'},'Interpreter','tex');
irf_legend(hca,{'R_{x}','R_{y}','R_{z}'},[0.88 0.10])

if 0
    hca = irf_panel('V_{MMS1}');
    irf_plot(hca,event.raw.Vi1_gsm);
    ylabel(hca,{'v_{GSM, MMS1}','(km/s)'},'Interpreter','tex');
    irf_legend(hca,{'v_{x}','v_{y}','v_{z}'},[0.88 0.10])
    
    mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
    hca = irf_panel('n_{p, MMS}'); set(hca,'ColorOrder',mmsColors)
    irf_pl_tx(hca,'event.raw.Ni?');
    ylabel(hca,{'n_p','(#/cm^{3})'},'Interpreter','tex');
    set(hca,'yscale','log'); % ,'log'
    irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.99 0.1],'color','cluster')
end
% current
hca = irf_panel('|J|'); % currently in micro Ampere
irf_plot(hca, event.raw.j_mag*1e6);
hold(hca,'on');
irf_plot(hca, event.nslv*1e6);
irf_plot(hca, event.Jmag_pks*1e6, 'linestyle', 'x');
hold(hca,'off');
grid on;
ylabel(hca, {'|J_{GSM}|','(\mu A/m^{2})'},'Interpreter','tex');
irf_legend(hca,{'data','noise level','peaks'},[0.88 0.90]);


% div over curl of B
if 1
    hca = irf_panel('divovercurl');
    irf_plot(hca,event.raw.divovercurl);
    ylabel(hca,{'|\nabla . B|/|\nabla \times B|'},'Interpreter','tex');
end

if 0
% J*r dirction
    hca = irf_panel('J*r all');
    irf_plot(hca, event.uv.jr_all, 'o', 'MarkerSize', 3);
    ylabel(hca, '$\hat{\mathbf{J}} \cdot \hat{\mathbf{r}}$', 'Interpreter', 'latex');
    
    hca = irf_panel('J*r pks');
    irf_plot(hca, event.uv.jr_pks, 'linestyle', 'o');
    ylabel(hca, '$(\hat{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}$', 'Interpreter', 'latex');
    
    % J*r ~direction and ~J_magnitude
    hca = irf_panel('J*r all vec');
    irf_plot(hca, event.jr_mag.all, 'o', 'MarkerSize', 3);
    ylabel(hca, '$\vec{\mathbf{J}} \cdot \hat{\mathbf{r}}$', 'Interpreter', 'latex');
    
    hca = irf_panel('J*r pks vec');
    irf_plot(hca, event.jr_mag.pks, 'linestyle', 'o');
    ylabel(hca, '$(\vec{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}$', 'Interpreter', 'latex');
end

irf_zoom(h(1:4),'x',tint_beob); % tint_beob [tint_beob(1)+300;tint_beob(1)+360]
irf_plot_axis_align(1,h(1:4));

tmp.filename = fullfile('events/', event.tint_string,'relevant_data');
irf_print_fig(tmp.filename,'png');
close(gcf);

%% same data for smaller tint
h = irf_figure(1,4,'reset');

hca = irf_panel('B_{MMS1}');
irf_plot(hca,event.raw.B1_gsm);
ylabel(hca,{'B_{GSM, MMS1}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10], 'fontsize', 14)

hca = irf_panel('R_{MMS1}');
irf_plot(hca,TSeries(event.raw.R1_gsm.time, event.raw.R1_gsm.data(:,1:3)/6378));
ylabel(hca,{'R_{GSM, MMS1}','(R_E)'},'Interpreter','tex');
irf_legend(hca,{'R_{x}','R_{y}','R_{z}'},[0.88 0.10], 'fontsize', 14)

% current
hca = irf_panel('|J|'); % currently in micro Ampere
irf_plot(hca, event.raw.j_mag*1e6);
hold(hca,'on');
irf_plot(hca, event.nslv*1e6);
irf_plot(hca, event.Jmag_pks*1e6, 'linestyle', 'x');
hold(hca,'off');
grid on;
ylabel(hca, {'|J_{GSM}|','(\mu A/m^{2})'},'Interpreter','tex');
irf_legend(hca,{'data';'noise level';'peaks'},[0.88 0.90], 'fontsize', 14);

% div over curl of B
if 1
    hca = irf_panel('divovercurl');
    irf_plot(hca,event.raw.divovercurl);
    ylabel(hca,{'|\nabla . B|/|\nabla \times B|'},'Interpreter','tex');
end

irf_plot_axis_align(1,h(1:4))
irf_zoom(h(1:4),'x',EpochTT([tint_beob(1)+300;tint_beob(1)+360]));


tmp.filename = fullfile('events/', event.tint_string,'relevant_data_zoom');
irf_print_fig(tmp.filename, 'png')
%irf_print_fig(tmp.filename,'png');
close(gcf);
%% spacecraft magnetic field measurements
h = irf_figure(1,3,'reset');

hca = irf_panel('B_{x}');
irf_plot(hca,event.raw.B1_gsm.x);
hold(hca,'on');
irf_plot(hca,event.raw.B2_gsm.x);
irf_plot(hca,event.raw.B3_gsm.x);
irf_plot(hca,event.raw.B4_gsm.x);
hold(hca,'off');
ylabel(hca,{'B_{x}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'MMS_{1}','MMS_{2}','MMS_{3}','MMS_{4}'},[0.88 0.10]);

hca = irf_panel('B_{y}');
irf_plot(hca,event.raw.B1_gsm.y);
hold(hca,'on');
irf_plot(hca,event.raw.B2_gsm.y);
irf_plot(hca,event.raw.B3_gsm.y);
irf_plot(hca,event.raw.B4_gsm.y);
hold(hca,'off');
ylabel(hca,{'B_{y}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'MMS_{1}','MMS_{2}','MMS_{3}','MMS_{4}'},[0.88 0.10]);


hca = irf_panel('B_{z}');
irf_plot(hca,event.raw.B1_gsm.z);
hold(hca,'on')
irf_plot(hca,event.raw.B2_gsm.z);
irf_plot(hca,event.raw.B3_gsm.z);
irf_plot(hca,event.raw.B4_gsm.z);
hold(hca,'off');
ylabel(hca,{'B_{z}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'MMS_{1}','MMS_{2}','MMS_{3}','MMS_{4}'},[0.88 0.10]);


irf_plot_axis_align(1,h(1:3))
irf_zoom(h(1:3),'x',[tint_beob(1)+300;tint_beob(1)+360]);

tmp.filename = fullfile('events/', event.tint_string,'B_each_sc');
irf_print_fig(tmp.filename,'png');
close(gcf);

%% spacecraft position plot
h = irf_figure(1,3,'reset');

hca = irf_panel('R_{x}');
irf_plot(hca,TSeries(event.raw.R1_gsm.time,event.raw.R1_gsm.data(:,1)./6378));
hold(hca,'on')
irf_plot(hca,TSeries(event.raw.R2_gsm.time,event.raw.R2_gsm.data(:,1)./6378));
irf_plot(hca,TSeries(event.raw.R3_gsm.time,event.raw.R3_gsm.data(:,1)./6378));
irf_plot(hca,TSeries(event.raw.R4_gsm.time,event.raw.R4_gsm.data(:,1)./6378));
hold(hca,'off');
ylabel(hca,{'R_{x}','$(R_{E})$'},'Interpreter','tex');
irf_legend(hca,{'MMS_{1}','MMS_{2}','MMS_{3}','MMS_{4}'},[0.88 0.10]);

hca = irf_panel('R_{y}');
irf_plot(hca,TSeries(event.raw.R1_gsm.time,event.raw.R1_gsm.data(:,2)./6378));
hold(hca,'on')
irf_plot(hca,TSeries(event.raw.R2_gsm.time,event.raw.R2_gsm.data(:,2)./6378));
irf_plot(hca,TSeries(event.raw.R3_gsm.time,event.raw.R3_gsm.data(:,2)./6378));
irf_plot(hca,TSeries(event.raw.R4_gsm.time,event.raw.R4_gsm.data(:,2)./6378));
hold(hca,'off');
ylabel(hca,{'R_{y}','$(R_{E})$'},'Interpreter','tex');
irf_legend(hca,{'MMS_{1}','MMS_{2}','MMS_{3}','MMS_{4}'},[0.88 0.10]);

hca = irf_panel('R_{z}');
irf_plot(hca,TSeries(event.raw.R1_gsm.time,event.raw.R1_gsm.data(:,3)./6378));
hold(hca,'on')
irf_plot(hca,TSeries(event.raw.R2_gsm.time,event.raw.R2_gsm.data(:,3)./6378));
irf_plot(hca,TSeries(event.raw.R3_gsm.time,event.raw.R3_gsm.data(:,3)./6378));
irf_plot(hca,TSeries(event.raw.R4_gsm.time,event.raw.R4_gsm.data(:,3)./6378));
hold(hca,'off');
ylabel(hca,{'R_{z}','$(R_{E})$'},'Interpreter','tex');
irf_legend(hca,{'MMS_{1}','MMS_{2}','MMS_{3}','MMS_{4}'},[0.88 0.10]);

irf_plot_axis_align(1,h(1:3))
irf_zoom(h(1:3),'x',[tint_beob(1)+300;tint_beob(1)+360]);

tmp.filename = fullfile('events/', event.tint_string,'R_each_sc');
irf_print_fig(tmp.filename,'png');
close(gcf);

%% current analyzing curlometer

h = irf_figure(1,5,'reset');

hca = irf_panel('B_{MMS1}');
irf_plot(hca,event.raw.B1_gsm);
ylabel(hca,{'B_{GSM, MMS1}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])

if 0
    hca = irf_panel('V_{MMS1}');
    irf_plot(hca,event.raw.Vi1_gsm);
    ylabel(hca,{'v_{GSM, MMS1}','(km/s)'},'Interpreter','tex');
    irf_legend(hca,{'v_{x}','v_{y}','v_{z}'},[0.88 0.10])
end
if 1
    mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
    hca = irf_panel('n_{p, MMS}'); set(hca,'ColorOrder',mmsColors)
    irf_pl_tx(hca,'event.raw.Ni?');
    ylabel(hca,{'n_p','(#/cm^{3})'},'Interpreter','tex');
    set(hca,'yscale','log'); % ,'log'
    irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.99 0.1],'color','cluster')
end
% current
hca = irf_panel('|J|'); % currently in micro Ampere
irf_plot(hca, event.raw.j_mag*1e6);
hold(hca,'on');
irf_plot(hca, event.nslv*1e6);
irf_plot(hca, event.Jmag_pks*1e6, 'linestyle', 'x');
hold(hca,'off');
grid on;
ylabel(hca, {'|J_{GSM}|','(\mu A/m^{2})'},'Interpreter','tex');
irf_legend(hca,{'data','noise level','peaks'},[0.88 0.90]);

if 0
    irf_plot(hca, event.Jmag_pks*1e6, 'linestyle', 'x');
    hold(hca,'off');
    grid on;
    ylabel(hca, {'|J_{GSM}|','(\mu A/m^{2})'},'Interpreter','tex');
    irf_legend(hca,{'data','noise level','peaks'},[0.88 0.90]);
end

% div over curl of B
if 0
    hca = irf_panel('divovercurl');
    irf_plot(hca,event.raw.divovercurl);
    ylabel(hca,{'|\nabla . B|/|\nabla \times B|'},'Interpreter','tex');
end

if 1
% J*r dirction
    hca = irf_panel('J*r all');
    irf_plot(hca, event.uv.jr_all, 'o', 'MarkerSize', 3);
    ylabel(hca, '$\hat{\mathbf{J}} \cdot \hat{\mathbf{r}}$', 'Interpreter', 'latex');
    
    hca = irf_panel('J*r pks');
    irf_plot(hca, event.uv.jr_pks, 'linestyle', 'o');
    ylabel(hca, '$(\hat{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}$', 'Interpreter', 'latex');
end
if 0
    % J*r ~direction and ~J_magnitude
    hca = irf_panel('J*r all vec');
    irf_plot(hca, event.jr_mag.all*1e6, 'o', 'MarkerSize', 3);
    ylabel(hca, {'\vec{\mathbf{J}} \cdot \hat{\mathbf{r}}','(\mu A/m^{2})'}, 'Interpreter', 'tex');
    
    hca = irf_panel('J*r pks vec');
    irf_plot(hca, event.jr_mag.pks*1e6, 'linestyle', 'o');
    ylabel(hca, {'\vec{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}','(\mu A/m^{2})'}, 'Interpreter', 'tex');
end

irf_plot_axis_align(1,h(1:5))
irf_zoom(h(1:5),'x',tint_beob); % tint_beob [tint_beob(1)+300;tint_beob(1)+360]

tmp.filename = fullfile('events/', event.tint_string,'Jr_uv');
irf_print_fig(tmp.filename,'png');
close(gcf);
%% current analyzing curlometer

h = irf_figure(1,5,'reset');

hca = irf_panel('B_{MMS1}');
irf_plot(hca,event.raw.B1_gsm);
ylabel(hca,{'B_{GSM, MMS1}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])

if 0
    hca = irf_panel('V_{MMS1}');
    irf_plot(hca,event.raw.Vi1_gsm);
    ylabel(hca,{'v_{GSM, MMS1}','(km/s)'},'Interpreter','tex');
    irf_legend(hca,{'v_{x}','v_{y}','v_{z}'},[0.88 0.10])
end
if 1
    mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
    hca = irf_panel('n_{p, MMS}'); set(hca,'ColorOrder',mmsColors)
    irf_pl_tx(hca,'event.raw.Ni?');
    ylabel(hca,{'n_p','(#/cm^{3})'},'Interpreter','tex');
    set(hca,'yscale','log'); % ,'log'
    irf_legend(hca,{'MMS1','MMS2','MMS3','MMS4'},[0.99 0.1],'color','cluster')
end
% current
hca = irf_panel('|J|'); % currently in micro Ampere
irf_plot(hca, event.raw.j_mag*1e6);
hold(hca,'on');
irf_plot(hca, event.nslv*1e6);
irf_plot(hca, event.Jmag_pks*1e6, 'linestyle', 'x');
hold(hca,'off');
grid on;
ylabel(hca, {'|J_{GSM}|','(\mu A/m^{2})'},'Interpreter','tex');
irf_legend(hca,{'data','noise level','peaks'},[0.88 0.90]);

if 0
    irf_plot(hca, event.Jmag_pks*1e6, 'linestyle', 'x');
    hold(hca,'off');
    grid on;
    ylabel(hca, {'|J_{GSM}|','(\mu A/m^{2})'},'Interpreter','tex');
    irf_legend(hca,{'data','noise level','peaks'},[0.88 0.90]);
end

% div over curl of B
if 0
    hca = irf_panel('divovercurl');
    irf_plot(hca,event.raw.divovercurl);
    ylabel(hca,{'|\nabla . B|/|\nabla \times B|'},'Interpreter','tex');
end

if 0
% J*r dirction
    hca = irf_panel('J*r all');
    irf_plot(hca, event.uv.jr_all, 'o', 'MarkerSize', 3);
    ylabel(hca, '$\hat{\mathbf{J}} \cdot \hat{\mathbf{r}}$', 'Interpreter', 'latex');
    
    hca = irf_panel('J*r pks');
    irf_plot(hca, event.uv.jr_pks, 'linestyle', 'o');
    ylabel(hca, '$(\hat{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}$', 'Interpreter', 'latex');
end
if 1
    % J*r ~direction and ~J_magnitude
    hca = irf_panel('J*r all vec');
    irf_plot(hca, event.jr_mag.all*1e6, 'o', 'MarkerSize', 3);
    ylabel(hca, {'$\vec{\mathbf{J}} \cdot \hat{\mathbf{r}}$','(\mu A/m^{2})'}, 'Interpreter', 'latex');
    
    hca = irf_panel('J*r pks vec');
    irf_plot(hca, event.jr_mag.pks*1e6, 'linestyle', 'o');
    ylabel(hca, {'$(\vec{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}$','(\mu A/m^{2})'}, 'Interpreter', 'latex');
end

irf_plot_axis_align(1,h(1:5))
irf_zoom(h(1:5),'x',tint_beob); % tint_beob [tint_beob(1)+300;tint_beob(1)+360]

tmp.filename = fullfile('events/', event.tint_string,'Jr_mag');
irf_print_fig(tmp.filename,'png');
close(gcf);
%% velocity distributions electrons

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
%%
%irf_cdf_read('*','*.cdf') % view raw data file to look at all fgm attributes
% B_DMPA
% B__GSE
% B_GSM
% B_BCS
% t
% R_GSE
% R_GSM
% errors and stupp
