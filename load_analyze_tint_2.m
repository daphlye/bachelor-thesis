function [omni, mmsd, event, data, current_info] = load_analyze_tint_2(tint, tint_beob, Param, i)
%LOAD_ANALYZE_TINT_2 downloads needed data and tests criteria
%   OUTPUT: omni data, mms data, only tint of 10min length, results of
%   testing current direction
%


%% load and analyze the tint 
    % load required OMNI and MMS data
    % load OMNI
    
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
    
    % load and resample (linear interpolation) MMS
    
    c_eval('mmsd.B?_gsm = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gsm_srvy_l2'',tint);',Param.SCall);
    % resampling on B since this has the highest resolution
    c_eval('mmsd.B?_gsm = mmsd.B?_gsm.resample(mmsd.B1_gsm);',2:4);
    c_eval('mmsd.R?_gsm = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_r_gsm_srvy_l2'',tint);', Param.SCall); % GSM [km] evaluates data for R 
    c_eval('mmsd.R?_gsm = mmsd.R?_gsm.resample(mmsd.R1_gsm);',2:4);
    c_eval('mmsd.Rres?  = mmsd.R?_gsm.resample(mmsd.B1_gsm);');
    % R mec (look at measuring rate)
    %c_eval('mmsd.R?_mec  = mms.db_get_ts(''mms?_mec_srvy_l2_epht89d'',''mms?_mec_r_gsm'',tint);', Param.SCall);
    % same resampling rate :(
    c_eval('mmsd.Vi?_gsm = irf_gse2gsm(mms.db_get_ts(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_bulkv_gse_fast'',    tint));',Param.SCall); % bulk velocity 
    c_eval('mmsd.Vi?_gsm = mmsd.Vi?_gsm.resample(mmsd.Vi1_gsm);',2:4);
    c_eval('mmsd.Vires?  = mmsd.Vi?_gsm.resample(mmsd.B1_gsm);');
    c_eval('mmsd.Ni?     = mms.db_get_ts(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_numberdensity_fast'',tint);',Param.SCall); % density [#]
    c_eval('mmsd.Ni?     = mmsd.Ni?.resample(mmsd.Ni1);',2:4);
    mmsd.ni = irf.ts_scalar(mmsd.Ni1.time,(mmsd.Ni1.data+mmsd.Ni2.data+mmsd.Ni3.data+mmsd.Ni4.data)/4);
    mmsd.ni = mmsd.ni.resample(mmsd.B1_gsm);
    
    % confirming mms is in tetrahedral formation
    [mmsd.tetr.volTensor,mmsd.tetr.R_Center,mmsd.tetr.dR1,mmsd.tetr.dR2,mmsd.tetr.dR3,mmsd.tetr.dR4,mmsd.tetr.L,mmsd.tetr.E,mmsd.tetr.P] = c_4_r(mmsd.R1_gsm.data(:,2:4),mmsd.R2_gsm.data(:,2:4),mmsd.R3_gsm.data(:,2:4),mmsd.R4_gsm.data(:,2:4)); %time is same as mmsd.R?_gsm
    tetr_E = mean(mmsd.tetr.E, 'omitnan'); tetr_P = mean(mmsd.tetr.P, 'omitnan');
    
    % Curlometer technique (Method from Schwartz et al. 1998)
    % J, divB, jxB, divTshear, divPb (from Example_MMS_B_E_J.m)
    if 1
        [mmsd.j,mmsd.divB,mmsd.B,mmsd.jxB,mmsd.divTshear,mmsd.divPb] = c_4_j('mmsd.Rres?','mmsd.B?_gsm');
        mmsd.divovercurl = mmsd.divB; % error on j in percentage of value
        mmsd.divovercurl.data = abs(mmsd.divovercurl.data)./mmsd.j.abs.data;
    end
    % Test if mean cone angle is quasi perpendicular (OMNI)
    cone_a = mean(acosd(omni.Bxyz.data(:,1)./sqrt(omni.Bxyz.data(:,1).^2+omni.Bxyz.data(:,2).^2+omni.Bxyz.data(:,3).^2)), 'omitnan');
    cone_a_ok = cone_a > 60 && cone_a < 120;
    % --> Maria told me to do this instead of looking at bow shock normal
    % direction, should be good enough criteria. 
    
    % load final tint data
    % extracting stable interval without first 15mins
    
    c_eval('event.raw.B?_gsm = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gsm_srvy_l2'',tint_beob);',Param.SCall);
    % resampling on B since this has the highest resolution
    c_eval('event.raw.B?_gsm = event.raw.B?_gsm.resample(mmsd.B1_gsm);',2:4);
    c_eval('event.raw.R?_gsm = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_r_gsm_srvy_l2'',tint_beob);', Param.SCall); % GSM [km] evaluates data for R 
    c_eval('event.raw.R?_gsm = event.raw.R?_gsm.resample(event.raw.R1_gsm);',2:4);
    c_eval('event.raw.Rres?  = event.raw.R?_gsm.resample(event.raw.B1_gsm);');
    % R mec (look at measuring rate)
    %c_eval('mmsd.R?_mec  = mms.db_get_ts(''mms?_mec_srvy_l2_epht89d'',''mms?_mec_r_gsm'',tint);', Param.SCall);
    % same resampling rate :(
    c_eval('event.raw.Vi?_gsm = irf_gse2gsm(mms.db_get_ts(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_bulkv_gse_fast'',    tint));',Param.SCall); % bulk velocity ions
    c_eval('event.raw.Vi?_gsm = event.raw.Vi?_gsm.resample(event.raw.Vi1_gsm);',2:4);
    c_eval('event.raw.Ve?_gsm = irf_gse2gsm(mms.db_get_ts(''mms?_fpi_fast_l2_des-moms'',''mms?_des_bulkv_gse_fast'',    tint));',Param.SCall); % bulk velocity electrons
    c_eval('event.raw.Ve?_gsm = event.raw.Vi?_gsm.resample(event.raw.Ve1_gsm);',2:4);
    c_eval('event.raw.Vires?  = event.raw.Vi?_gsm.resample(event.raw.B1_gsm);');
    c_eval('event.raw.Ni?     = mms.db_get_ts(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_numberdensity_fast'',tint);',Param.SCall); % density [#]
    c_eval('event.raw.Ni?     = event.raw.Ni?.resample(event.raw.Ni1);',2:4);
    event.raw.ni = irf.ts_scalar(event.raw.Ni1.time,(event.raw.Ni1.data+event.raw.Ni2.data+event.raw.Ni3.data+event.raw.Ni4.data)/4);
    event.raw.ni = event.raw.ni.resample(event.raw.B1_gsm);
    
    %% Load PDist using mms.make_pdist
    if 0
        % des = dual electron spectrometer
        tmp.filepath_and_filename = mms.get_filepath(irf_ssub('mms?_fpi_fast_l2_des-dist',Param.ic),tint_beob);
        c_eval('[event.ePDist?,event.ePDistError?] = mms.make_pdist(tmp.filepath_and_filename);', Param.ic)
        % dis = dual ion spectrometer
        tmp.filepath_and_filename = mms.get_filepath(irf_ssub('mms?_fpi_fast_l2_dis-dist',Param.ic),tint_beob);
        c_eval('[event.iPDist?,event.iPDistError?] = mms.make_pdist(tmp.filepath_and_filename);', Param.ic)
            % load supporting data
        %   SC potential
        c_eval('event.scPot?=mms.db_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_scpot_fast_l2'',tint);', Param.ic);
        %   MMS EDP Fast L2 DCE (electric field), DSL (despun spacecraft level)
        c_eval('event.dslE?=mms.db_get_ts(''mms?_edp_fast_l2_dce'',''mms?_edp_dce_dsl_fast_l2'',tint);', Param.ic);  
        
    end

    %% Curlometer technique (Method from Schwartz et al. 1998)
    % J, divB, jxB, divTshear, divPb (from Example_MMS_B_E_J.m)
    
    [event.raw.f.j,event.raw.f.divB,event.raw.f.B,event.raw.f.jxB,event.raw.f.divTshear,event.raw.f.divPb] = c_4_j('event.raw.Rres?','event.raw.B?_gsm');
    % uses first three rows for x y z. (confirmed looking at c_4_j.m >
    % c_4_grad.m)

    % extract curlometer datapoints where divovercurl is < 1
    event.raw.divovercurl = event.raw.f.divB;     event.raw.divovercurl.data = abs(event.raw.divovercurl.data)./event.raw.f.j.abs.data; % error on j in percentage of value

    tmp.ind = find(event.raw.divovercurl.data < 1);
    event.raw.j=TSeries(event.raw.f.j.time(tmp.ind),event.raw.f.j.data(tmp.ind,:));
    event.raw.divB=TSeries(event.raw.f.divB.time(tmp.ind),event.raw.f.divB.data(tmp.ind));
    event.raw.B=TSeries(event.raw.f.B.time(tmp.ind),event.raw.f.B.data(tmp.ind));
    event.raw.jxB=TSeries(event.raw.f.jxB.time(tmp.ind),event.raw.f.jxB.data(tmp.ind));
    event.raw.divTshear=TSeries(event.raw.f.divTshear.time(tmp.ind),event.raw.f.divTshear.data(tmp.ind));
    event.raw.divPb=TSeries(event.raw.f.divPb.time(tmp.ind),event.raw.f.divPb.data(tmp.ind));
    %   change position array to size of j array
    c_eval('event.raw.Rres?  = event.raw.R?_gsm.resample(event.raw.j);');
    
%% filter out "noise"
    event.raw.j_mag.data = sqrt(event.raw.j.data(:,1).^2 + event.raw.j.data(:,2).^2 + event.raw.j.data(:,3).^2);
    event.raw.j_mag = TSeries(event.raw.j.time, sqrt(event.raw.j.data(:,1).^2 + event.raw.j.data(:,2).^2 + event.raw.j.data(:,3).^2));
    [event.avgj, event.medj, event.stdj] = get_stats(event.raw.j_mag.data);
    % Find peaks in Jmag and create corresponding timeseries
    nslv = event.medj + 2*event.stdj; 
    %nslv = 2*event.stdj; 
    pkprom = 0.07*1e-6;
    [event.pksj,event.locsj]=findpeaks(event.raw.j_mag.data,'MinPeakHeight',event.stdj*2, 'MinPeakProminence',pkprom); %'MinPeakProminence',0.07*1e-6
    % building TSerieses with only j peaks data
    % j_mag
    tmp.D_pks = zeros(size(event.raw.j_mag.data));
    tmp.D_pks(:) = nan; tmp.D_pks(event.locsj) = event.pksj; 
    
    event.Jmag_pks = TSeries(event.raw.j_mag.time, tmp.D_pks); %j magnitude peaks
    event.j_pks = TSeries(event.raw.j.time(event.locsj), event.raw.j.data(event.locsj,:)); %j peaks
    
    % noise level definition from salinas:
    event.nslv = TSeries(event.raw.j_mag.time, (2*event.stdj)*ones(length(event.raw.j_mag.time), 1)); % noise level definition from salinas
    
    %% Particle moments method
    % calculate J = n_e*e*(v_i-v_e)
    event.J_pmm = event.raw.Ni1*1e6*1.602*1e-19*(event.raw.Vi1_gsm - event.raw.Ve1_gsm)*1e3;
    e = 1.602*10^(-19);
    %event.d_J_pmm = ((e*(event.raw.Vi1_gsm - event.raw.Ve1_gsm)*d_n_p)^2+(event.raw.Ni1*e*d_v_p)^2+(event.raw.Ni1*e*d_v_e)^2)^0.5;
    %% J*r stuff
    % calculating J*r direction
    %calculating unit vector of j and r     confirm where fgm r comes from!!! otherwise use mec
    event.uv.j_all = TSeries(event.raw.j.time, unit_vector(event.raw.j.data));
    event.uv.j_pks = TSeries(event.j_pks.time, unit_vector(event.j_pks.data));
    % calculating center of tetrahedra spacecraft from resampled R data
    [event.tetr.volTensor,event.tetr.R_Center, event.tetr.dR1,event.tetr.dR2,event.tetr.dR3,event.tetr.dR4,event.tetr.L,event.tetr.E,event.tetr.P] = c_4_r(event.raw.Rres1.data(:,2:4),event.raw.Rres2.data(:,2:4),event.raw.Rres3.data(:,2:4),event.raw.Rres4.data(:,2:4));
    event.uv.r = TSeries(event.raw.Rres1.time, unit_vector(event.tetr.R_Center)); 
    
    % calculating J*r
    %   case  1 - parallel, current towards earth (since neg Jx means
    %   towards earth - it is counterintuitive)
    %   case  0 - orthogonal current in ms directions
    %   case -1 - antiparallel current away from earth 

    % right direction
    if mean(event.raw.R1_gsm.data(:,2))> 0 % dawnside +y, J*r should point in +r, away from earth
        rd = 1;
    else
        rd = -1;
    end
    
    % without j noiselevelstuff (including all j data)
    %   dot product for that data set: elemetwise multiplication, rowwise summation
    event.uv.jr_all = TSeries(event.uv.j_all.time, sum(event.uv.j_all.data .* event.uv.r.data, 2));
    % with noiselevel limitation (only peaks)
    event.uv.jr_pks = TSeries(event.uv.j_pks.time, sum(event.uv.j_pks.data .* event.uv.r.data(event.locsj), 2));
    % calculating J*r direction with dependance on J-magnitude
    event.jr_mag.all = TSeries(event.raw.j.time, sum(event.raw.j.data .* event.uv.r.data, 2));
    % with noiselevel limitation (only peaks)
    event.jr_mag.pks = TSeries(event.j_pks.time, sum(event.j_pks.data .* event.uv.r.data(event.locsj), 2));
    %% analyzing
    position = [mean(event.raw.R1_gsm.data(:,1)./6378), mean(event.raw.R1_gsm.data(:,2)./6378), mean(event.raw.R1_gsm.data(:,3)./6378)]; 
    % J*r direction
    j_uv_all = jr_data_ana(event.uv.jr_all.data, event.raw.R1_gsm.data);
    j_uv_pks = jr_data_ana(event.uv.jr_pks.data, event.raw.R1_gsm.data);

    % Integrate over all j_maguv to see in which direction the overall
    % current goes
    tot_curr_all = sum(event.jr_mag.all.data); % Sum of all current values
    tot_curr_pks = sum(event.jr_mag.pks.data);

    % summing up over magnitude of currents that are pointing towards / away (above
    % |0.25|) 
    tot_curr_a25_all = sum(event.jr_mag.all.data(find(event.uv.jr_all.data > 0.25)))+ sum(event.jr_mag.all.data(find(event.uv.jr_all.data < -0.25)));
    tot_curr_a25_pks = sum(event.jr_mag.pks.data(find(event.uv.jr_pks.data > 0.25)))+sum(event.jr_mag.pks.data(find(event.uv.jr_pks.data < -0.25)));

    %% putting everything in a table
    data = [cone_a, nslv, position, tetr_E, tetr_P];
    current_info = [rd, j_uv_all(:,6),j_uv_pks(:,6), tot_curr_all*1e6, tot_curr_pks*1e9, tot_curr_a25_all*1e6, tot_curr_a25_pks*1e9];
end