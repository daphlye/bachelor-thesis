function [plt] = plot_tint_3(omni, mmsd, event, tint, tint_beob, Param)
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
% Set default legend properties
set(groot, 'defaultLegendBox', 'off');               % Remove box
set(groot, 'defaultLegendColor', 'none');            % Transparent background
set(groot, 'defaultLegendFontSize', 10);             % Smaller font size
set(groot, 'defaultLegendLocation', 'best');         % Automatically position

%plt = 0; if time, we can make this a finction that tests if all plots
%got created so that code can stop if not.
%% test
plot_omni_mms;

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
legend(hca,{'B_{x}','B_{y}','B_{z}'},'location', 'best')

if 0 
hca = irf_panel('R_{MMS1}');
    irf_plot(hca,TSeries(event.raw.R1_gsm.time, (event.raw.R1_gsm.data(:,1:3)/6378)));
    ylabel(hca,{'R_{GSM, MMS1}','(R_E)'},'Interpreter','tex');
    legend(hca,{'R_{x}','R_{y}','R_{z}'},'location')
end

if 1
    hca = irf_panel('V_{MMS1}');
    irf_plot(hca,event.raw.Vi1_gsm);
    ylabel(hca,{'v_{GSM, MMS1}','(km/s)'},'Interpreter','tex');
    irf_legend(hca,{'v_{x}','v_{y}','v_{z}'},[0.88 0.10])
    
    mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
    hca = irf_panel('n_{p}'); set(hca,'ColorOrder',mmsColors)
    irf_plot(hca,event.raw.Ni1);
    ylabel(hca,{'n_{p}','(#/cm^{3})'},'Interpreter','tex');
end
% current
hca = irf_panel('|J|'); % currently in micro Ampere
irf_plot(hca, event.raw.j_mag);
hold(hca,'on');
irf_plot(hca, event.nslv);
irf_plot(hca, event.Jmag_pks, 'linestyle', 'x');
hold(hca,'off');
grid on;
ylabel(hca, {'|J_{GSM}|','(\mu A/m^{2})'},'Interpreter','tex');
legend(hca,{'data','noise level','peaks'},'location', 'best');


% div over curl of B
if 1
    hca = irf_panel('divovercurl');
    irf_plot(hca,event.raw.divovercurl);
    ylabel(hca, '$\frac{|\nabla \cdot \mathbf{B}|}{|\nabla \times \mathbf{B}|}$', 'interpreter', 'latex'); % 
end

linkaxes(h,'x');
irf_zoom(h(1:4),'x',tint_beob); % tint_beob [tint_beob(1)+300;tint_beob(1)+360]

tmp.filename = fullfile('events/plots/','relevant_data/');
if ~exist(tmp.filename, 'dir')
    mkdir(tmp.filename);
end
tmp.filename = fullfile('events/plots/relevant_data/', event.tint_string);
irf_print_fig(tmp.filename,'png');
close(gcf);

%% current analyzing curlometer

h = irf_figure(1,5,'reset');

hca = irf_panel('B_{MMS1}');
irf_plot(hca,event.raw.B1_gsm);
ylabel(hca,{'B_{GSM}','(nT)'},'interpreter','tex');
legend(hca, {'$B_{x}$', '$B_{y}$', '$B_{z}$'}, 'interpreter', 'latex');

if 0
    hca = irf_panel('V_{MMS1}');
    irf_plot(hca,event.raw.Vi1_gsm);
    ylabel(hca,{'v_{GSM}','(km/s)'},'Interpreter','tex');
    irf_legend(hca,{'v_{x}','v_{y}','v_{z}'},[0.88 0.10])
end
if 1
    mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
    hca = irf_panel('n_{p}'); set(hca,'ColorOrder',mmsColors)
    irf_plot(hca,event.raw.Ni1);
    ylabel(hca,{'n_{p}','(#/cm^{3})'},'Interpreter','tex');
    %set(hca,'yscale','log'); % ,'log'
    %legend(hca,{'MMS1','MMS2','MMS3','MMS4'},'Location', 'best')
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
legend(hca,{'data','noise level','peaks'},'location','north');

if 1
% J*r dirction
    hca = irf_panel('J*r all');
    irf_plot(hca, event.uv.jr_all, 'o', 'MarkerSize', 3);
    ylabel(hca, '$\hat{\mathbf{J}} \cdot \hat{\mathbf{r}}$', 'interpreter', 'latex');
    
    hca = irf_panel('J*r pks');
    irf_plot(hca, event.uv.jr_pks, 'linestyle', 'o');
    ylabel(hca, '$(\hat{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}$', 'interpreter', 'latex');
end

linkaxes(h,'x');
%irf_plot_axis_align(1,h(1:5))
irf_zoom(h(1:5),'x',tint_beob); % tint_beob [tint_beob(1)+300;tint_beob(1)+360]

tmp.filename = fullfile('events/plots/','Jr_uv/');
if ~exist(tmp.filename, 'dir')
    mkdir(tmp.filename);
end
tmp.filename = fullfile('events/plots/','Jr_uv/', event.tint_string);
irf_print_fig(tmp.filename,'png');
close(gcf);
%% current analyzing curlometer

h = irf_figure(1,5,'reset');

hca = irf_panel('B_{MMS1}');
irf_plot(hca,event.raw.B1_gsm);
ylabel(hca,{'B_{GSM, MMS1}','(nT)'},'interpreter','tex');
legend(hca,{'B_{x}','B_{y}','B_{z}'},'Location', 'best')

if 1
    mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
    hca = irf_panel('n_{p}'); set(hca,'ColorOrder',mmsColors)
    irf_plot(hca,event.raw.Ni1);
    ylabel(hca,{'n_{p}','(#/cm^{3})'},'Interpreter','tex');
end
% current
hca = irf_panel('|J|'); % currently in micro Ampere
irf_plot(hca, event.raw.j_mag*1e6);
hold(hca,'on');
irf_plot(hca, event.nslv*1e6);
irf_plot(hca, event.Jmag_pks*1e6, 'linestyle', 'x');
hold(hca,'off');
grid on;
ylabel(hca, {'|J_{GSM}|','(\mu A/m^{2})'},'interpreter','tex');
legend(hca,{'data','noise level','peaks'},'Location', 'best');

if 1
    % J*r ~direction and ~J_magnitude
    hca = irf_panel('J*r all vec');
    irf_plot(hca, event.jr_mag.all*1e6, 'o', 'MarkerSize', 3);
    ylabel(hca, {'$\vec{\mathbf{J}} \cdot \hat{\mathbf{r}}$','($\mu A/m^{2}$)'}, 'interpreter', 'latex');
    
    hca = irf_panel('J*r pks vec');
    irf_plot(hca, event.jr_mag.pks*1e6, 'linestyle', 'o');
    ylabel(hca, {'$(\vec{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}$','($\mu A/m^{2}$)'}, 'interpreter', 'latex');
end

irf_plot_axis_align(1,h(1:5))
irf_zoom(h(1:5),'x',tint_beob); % tint_beob [tint_beob(1)+300;tint_beob(1)+360]

tmp.filename = fullfile('events/plots/','Jr_mag/');
if ~exist(tmp.filename, 'dir')
    mkdir(tmp.filename);
end
tmp.filename = fullfile('events/plots/','Jr_mag/', event.tint_string);
irf_print_fig(tmp.filename,'png');
close(gcf);


