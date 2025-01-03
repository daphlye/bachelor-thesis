% make tables 
tmp.filename = fullfile('events/tables');
if ~exist(tmp.filename, 'dir')
    mkdir(tmp.filename);
end

% Put information from arrays into readable time format
tmp.adjusted_time = good_time(:,1) + (15 * 60);
tmp.adjusted_time = datetime(tmp.adjusted_time, 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
tmp.end_time = datetime(good_time(:,2), 'ConvertFrom', 'posixtime', 'TimeZone', 'UTC');
values.t_start = tmp.adjusted_time; % Assign datetime array
values.t_end = tmp.end_time;        % Assign datetime array
values{:,3:end} = val_arr;
current{:, :} = curr_arr;
values_time = array2table(zeros(height(good_time),2));
values_time.Properties.VariableNames = {'t_start', 't_end'};
values_time = array2table(values{:,1:2});
valtime_ltx = table2latex(values_time);

%change j?noiselevel to muA/m^2
values{:,4} = values{:,4}*1e6;

values_ltx = table2latex(values);
fid = fopen('events/tables/values_ltx.txt', 'w');  % Open the file for writing
fprintf(fid, '%s\n', values_ltx);           % Write each line of text
fclose(fid);  
current_ltx = array2latex(current);
fid = fopen('events/tables/current_ltx.txt', 'w');  % Open the file for writing
fprintf(fid, '%s\n', current_ltx);           % Write each line of text
fclose(fid);  
errors_ltx = array2latex(errors_for_table);

% list event numbers tints
    event_tints = array2table(zeros(height(good_time),3));
    event_tints.Properties.VariableNames = {'Event','Start','End'};
    event_tints{:,1}= (1:length(good_time))';
    event_tints{:, 2} = string(values{:,1});
    event_tints{:, 3} = string(values{:,2});
    % make table
    times_ltx = table2latex(event_tints);
    fid = fopen('events/tables/event_tint.txt', 'w');  % Open the file for writing
    fprintf(fid, '%s\n', times_ltx);           % Write each line of text
    fclose(fid);  


%% calculate Jpks_>0.25/Jpks
parts_ab25 = curr_arr(:,7)./curr_arr(:,5);
newcurr_arr = zeros(length(curr_arr), 6);
newcurr_arr(:,1:5) = curr_arr(:,1:5);
newcurr_arr(:,6) = parts_ab25*100;
% calculate error
Jpks25 = curr_arr(:,7); 
Jpks   = curr_arr(:,5);
errJpks25 = errors_for_table(:,3);
errJpks   = errors_for_table(:,2);

error.Jab25perJ = sqrt((errJpks25./Jpks).^2 + (Jpks25./(Jpks.^2).*errJpks).^2);
error.newcurr_arr = [zeros(10,3), errors_for_table(:,1:2),error.Jab25perJ];
% calculate mean of arrays for all values and the related error
% Initialize arrays to store column means and propagated errors
colMeans = zeros(1, size(newcurr_arr, 2));
colErrors = zeros(1, size(error.newcurr_arr, 2));

% Iterate through each column
for col = 1:size(newcurr_arr, 2)
    % Extract values and errors for the current column
    colData = newcurr_arr(:, col);
    colError = error.newcurr_arr(:, col);
    
    % Calculate the mean of the column
    colMeans(col) = mean(colData);
    
    % Propagate the error for the mean
    colErrors(col) = sqrt(sum(colError.^2)) / (length(colData).^2);
end

%% create matrix for color 
comat_arr = newcurr_arr(:,2:end);      
minusy = find(newcurr_arr(:,1)==-1);
comat_arr(minusy,3:4) = newcurr_arr(minusy,4:5)*(-1); % changing sign of values on -y side for unified colormap

colorMatrix = array2table(comat_arr);
eventMatrix = newcurr_arr(:,2:end);   % Matrix for displaying numbers (with correct sign of y)
criteriaLabels = {'(a) $(\mathbf{\hat{J}}\cdot\mathbf{\hat{r}})_{all} [\%]$', '(b) $(\mathbf{\hat{J}}\cdot\mathbf{\hat{r}})_{peaks} [\%]$', '(c) $(\mathbf{\vec{J}}\cdot\mathbf{\hat{r}})_{all} [\mu A/ m^2]$', ...
    '(d) $(\mathbf{\vec{J}}\cdot\mathbf{\hat{r}})_{peaks} [nA/ m^2]$', '(e) $(\frac{\mathbf{\vec{J}}\cdot\mathbf{\hat{r}}_{> 0.25}}{\mathbf{\vec{J}}\cdot\mathbf{\hat{r}}})_{peaks} [\%]$'};
eventLabels = 1:length(good_time);

%manually exclude event 9 because current values exceed 
comatwo13 = comat_arr;
comatwo13(9, :) = 0;

% create colorbar limits for each column
limit_perc = [0,100];
limit_all  = [max(abs(comatwo13(:,3)))*-1,max(abs(comatwo13(:,3)))];
limit_pks  = [max(abs(comatwo13(:,4)))*-1,max(abs(comatwo13(:,4)))];
limit_perc_special = [-100,100]; % includes negative percentage because 
% currents above 0.25 might point in different direction than mean

colorlimitsArray = {limit_perc, limit_perc, limit_all, limit_pks, limit_perc_special}; % limit_perc, limit_perc, limit_all, limit_pks, limit_alla25, limit_pksa25
% make correlation table
correlation_matrix_pmy(colorMatrix, eventMatrix, eventLabels, criteriaLabels, colorlimitsArray)


%% calculating mean current density 
mean_curr = mean(newcurr_arr(:, 5));

%% showcase interval of good SW situation 
tint2 = good_time(2, :); %irf.tint('2024-12-12T08:00:00.000Z/2024-12-12T10:00:00.000Z');
tint = [tint2(1)-60*60,tint2(1)+60*90]; % take 2 hour tint around first interval of many or so. 

    tmp.B = irf_get_data_omni(tint, 'b,bx,byGSM,bzGSM', 'omni_min');% loads 1-min OMNI data flow speed velocity in gsm coordinates (OMNI v data is in gse) for selected tint
    tmp.V = irf_gse2gsm(irf_get_data_omni(tint, 'v,vx,vy,vz', 'omni_min'));
    tmp.Ms_n = irf_get_data_omni(tint, 'Ms,n', 'omni_min'); % 1 AU IP Magnetosonic Mach number, ion dens
    tmp.IND = find(tmp.Ms_n(:,2)>99);    tmp.Ms_n(tmp.IND,2) = NaN; % removing fillvalues
    omni.time = EpochUnix(tmp.B(:,1:1)); % convert downloaded time to unix data
    omni.B_tot = TSeries(omni.time, [tmp.B(:,2:2)], 'to', 1); % GSE [nT]
    omni.Bxyz = TSeries(omni.time, [tmp.B(:,3:5)], 'to', 1); % GSM [nT]
    omni.V_tot = TSeries(omni.time, [tmp.V(:,2:2)], 'to', 1); % GSM [km/s]
    omni.Vxyz = TSeries(omni.time, [tmp.V(:,3:5)], 'to', 1); % GSM [km/s]
    omni.Ms = TSeries(omni.time, [tmp.Ms_n(:,2:2)], 'to', 1); % [#]
    omni.np = TSeries(omni.time, [tmp.Ms_n(:,3:3)], 'to', 1); % [cc]
%% calculating means and maxes that are needed 
meanomni.sw.B = [mean(omni.Bxyz.data(:,1)), mean(omni.Bxyz.data(:,2)), mean(omni.Bxyz.data(:,3)), mean(omni.B_tot.data)];
meanomni.sw.B_std = [std(omni.Bxyz.data(:,1)), std(omni.Bxyz.data(:,2)), std(omni.Bxyz.data(:,3)), std(omni.B_tot.data)];

meanomni.sw.V = [mean(omni.Vxyz.data(:,1),'omitnan'), mean(omni.Vxyz.data(:,2),'omitnan'), mean(omni.Vxyz.data(:,3),'omitnan'), mean(omni.V_tot.data,'omitnan')];
meanomni.sw.V_std = [std(omni.Vxyz.data(:,1),'omitnan'), std(omni.Vxyz.data(:,2),'omitnan'), std(omni.Vxyz.data(:,3),'omitnan'), std(omni.V_tot.data,'omitnan')];
meanomni.sw.N = mean(omni.np.data,'omitnan');
meanomni.sw.Ms = mean(omni.Ms.data,'omitnan');
maxcurr = max(event.raw.j_mag.data)*1E9;
maxcurr_err = error.J(find(max(event.raw.j_mag.data)))*1E9;

%%
h = irf_plot(4,'newfigure');

hca = irf_panel('M_{MS, OMNI}');
irf_plot(hca,omni.Ms);
ylabel(hca,{'M_{MS}','(#)'},'Interpreter','tex');
irf_zoom(hca,'x',tint);
irf_zoom(hca,'y');
%irf_legend(hca,'(a)',[0.99 0.98],'color','k')

hca = irf_panel('B_{OMNI}');
irf_plot(hca,omni.Bxyz);
ylabel(hca,{'B_{GSM}','(nT)'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
%irf_zoom(hca,'x',tint);
%irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca = irf_panel('v_{OMNI}');
irf_plot(hca,omni.Vxyz);
ylabel(hca,{'v_{GSM}','(km/s)'},'Interpreter','tex');
irf_legend(hca,{'v_{x}','v_{y}','v_{z}'},[0.88 0.10])
%irf_zoom(hca,'x',tint);
%irf_legend(hca,'(c)',[0.99 0.98],'color','k')

hca = irf_panel('ni_{OMNI}'); % set(hca,'ColorOrder',mmsColors)
irf_plot(hca,omni.np);
ylabel(hca,{'n_{p}','(#/cm^{3})'},'Interpreter','tex'); % Units confirmed (SPDF)
%set(hca,'yscale','log');
%irf_zoom(hca,'y',[0 10]);
%irf_legend(hca,'(d)',[0.99 0.98],'color','k')
irf_plot_axis_align(1,h(1:4))
irf_zoom(h(1:4),'x',tint);
irf_pl_mark(h,tint2)
irf_pl_number_subplots
%saveas(h, 'plots/tint1/omni');  % Save the current figure as a .fig file
tmp.filename = fullfile('events/plots/','omni/');
if ~exist(tmp.filename, 'dir')
    mkdir(tmp.filename);
end
tmp.filename = fullfile('events/plots/omni/',event.tint_string);
irf_print_fig(tmp.filename,'png');
%%
h = irf_figure(1,5,'reset');
% B components
hca = irf_panel('B_{MMS1}');
irf_plot(hca,event.raw.B1_gsm);
ylabel(hca,{'B_{x,y,z, GSM}','(nT)'},'interpreter','tex');
legend(hca,{'B_{x}','B_{y}','B_{z}'},'Location', 'southeast')
% B magnitude
magnitude = hypot(hypot(event.raw.B1_gsm.data(:,1), event.raw.B1_gsm.data(:,2)), event.raw.B1_gsm.data(:,3));
B_mag = TSeries(event.raw.B1_gsm.time,magnitude);
hca = irf_panel('B_{mag}');
irf_plot(hca,B_mag);
ylabel(hca,{'B_{tot, GSM}','(nT)'},'interpreter','tex');
% n
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
legend(hca,{'data','noise level','peaks'},'Location', 'east');

if 1
    hca = irf_panel('J*r pks vec');
    irf_plot(hca, event.jr_mag.pks*1e6, 'linestyle', 'o');
    ylabel(hca, {'$(\vec{\mathbf{J}} \cdot \hat{\mathbf{r}})_{peaks}$','($\mu A/m^{2}$)'}, 'interpreter', 'latex');
   
end

irf_plot_axis_align(1,h(1:5))
tint3 = irf.tint('2015-12-31T23:28:00.000Z/2015-12-31T23:34:00.000Z');
irf_zoom(h(1:5),'x',tint3);
irf_pl_mark(h,tint_exp)
tint_exp = irf.tint('2015-12-31T23:30:00.000Z/2015-12-31T23:32:00.000Z');
irf_pl_mark(h,tint_exp)
irf_pl_number_subplots
