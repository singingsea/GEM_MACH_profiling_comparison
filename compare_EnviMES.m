function compare_EnviMES(target,instrument,instrument_name,comparison_with)
% this function is to compare EnviMes results with Pandora and AOEROCAN
save_fig = 1;
plot_path = 'D:\Projects\EnviMes\plot\';mkdir(plot_path);
DU = 2.6870e+16;
addpath(genpath('C:\Users\ZhaoX\Documents\MATLAB\matlab'));

% target = 'NO2_428NM';instrument = '1694_2';instrument_name = 'EnviMes';comparison_with = 'DS_trop_column';%
% target = 'NO2_460NM';instrument = '1694_1';instrument_name = 'EnviMes';comparison_with = 'DS_trop_column';%
% target = 'NO2_428NM';instrument = '1694_2';instrument_name = 'EnviMes';comparison_with = 'PGN_trop_column';%
% target = 'NO2_460NM';instrument = '1694_1';instrument_name = 'EnviMes';comparison_with = 'PGN_trop_column';%
% target = 'NO2_428NM';instrument = '1694_2';instrument_name = 'EnviMes';comparison_with = 'PGN_surface';%
% target = 'NO2_460NM';instrument = '1694_1';instrument_name = 'EnviMes';comparison_with = 'PGN_surface';%
% target = 'NO2_428NM';instrument = '1694_2';instrument_name = 'EnviMes';comparison_with = 'in_situ';%
% target = 'NO2_460NM';instrument = '1694_1';instrument_name = 'EnviMes';comparison_with = 'in_situ';%
% target = 'O4_360NM';instrument = '1694_2';instrument_name = 'EnviMes';comparison_with = 'AEROCAN';%
% target = 'O4_477NM';instrument = '1694_1';instrument_name = 'EnviMes';comparison_with = 'AEROCAN';%

% target = 'O4_360NM';instrument = '1694_2';instrument_name = 'EnviMes';comparison_with = 'in_situ_PM';%
% target = 'O4_477NM';instrument = '1694_1';instrument_name = 'EnviMes';comparison_with = 'in_situ_PM';%
% leinnterp1';% AEROCAN interp to 360 nm
% target = 'O4_477NM';instrument = '1694_1';instrument_name = 'EnviMes';comparison_with = 'AEROCAN interp2';% AEROCAN interp to 477 nm

% target = 'HCHO_343NM';instrument = '1694_2';instrument_name = 'EnviMes';comparison_with = 'FTIR';% 
% target = 'HCHO_343NM';instrument = '1694_2';instrument_name = 'EnviMes';comparison_with = 'FTIR_surface';% 



plot_path = ['D:\Projects\EnviMes\plot\' comparison_with '\'];mkdir(plot_path);


% target = 'NO2_428NM';instrument = '1694_2';instrument_name = 'EnviMes';% tragets = {'NO2_428NM','O4_360NM','HCHO_343NM'};
% target = 'O4_477NM';instrument = '1694_1';instrument_name = 'EnviMes';% tragets = {'NO2_460NM','O4_477NM'};
% comparison_with = 'PGN_surface';
% comparison_with = 'in_situ';
% comparison_with = 'AEROCAN';
%% sort out comparision data location

EnviMes_data_path = ['D:\Projects\EnviMes\output_profiles\' instrument '\'];
        
list = dir([EnviMes_data_path target '*']);
Env = [];
for i = 1:numel(list)
    disp(['Loading: ' list(i).name]);
    load([EnviMes_data_path list(i).name]);
    Env = [Env;data];
end
tf = Env.qa>= 2;% simple QC 
Env(tf,:) = [];
if strcmp(comparison_with, 'PGN_surface') | strcmp(comparison_with, 'in_situ') | strcmp(comparison_with, 'FTIR_surface')
    Env.vmr = Env.profile./Env.air_in_column.*200*100.*1e9;% convert from molec/cm^3 to ppbv; here the vertical grid height of MMF is 200 m
end

use_yyaxis = false;% default is use same y axis for two instrument/dataset, but for some comparison, e.g., PM2.5 vs. extincion, we need two y axis
switch comparison_with
    case 'DS_trop_column'
        load('C:\Projects\Pandora\output\103\PanPS_BlickP_mearged_NO2\estimated_Trop_plots\P103_DS_NO2_hq_offsetremoved_trop_bias_removed.mat');
        Second_dataset = PandoraNO2;      
        Second_instrument_name = 'Pandora';
        binwidth = 0.1;
        Second_dataset.target2 = Second_dataset.VCD_trop;% the PGN DS total column subtract OMI strat
        Env.target = Env.tvcd./DU;% EnviMes trop NO2
        measurand_txt = 'NO_2_t_r_o_p [DU]';% absolute observation
        measurand_txt2 = '(NO_2_t_r_o_p ratio)';% ratio text
        measurand_txt_target2 = measurand_txt;
        avg_time = 10; % [min] EnviMes and Pandora data will be averged into this time bin
        yrange = [0 1.2];
    case 'PGN_trop_column'
        load('C:\Projects\Pandora\output\103\PGN\Pandora103s1_Downsview_L2Trop_rnvh1p1-7_hq.mat');% Pandora multi-axis results from PGN
        Second_dataset = Pandora;      
        Second_instrument_name = 'Pandora';
        binwidth = 0.1;
        Second_dataset.target2 = Second_dataset.NO2_VCD_trop;% the PGN multi-axis data product
        Env.target = Env.tvcd./DU;% EnviMes trop NO2
        measurand_txt = 'NO_2_t_r_o_p [DU]';% absolute observation
        measurand_txt2 = '(NO_2_t_r_o_p ratio)';% ratio text
        measurand_txt_target2 = measurand_txt;
        avg_time = 10; % [min] EnviMes and Pandora data will be averged into this time bin
        yrange = [0 1.2];
    case 'PGN_surface'
        load('C:\Projects\Pandora\output\103\PGN\Pandora103s1_Downsview_L2Trop_rnvh1p1-7_hq.mat');% Pandora multi-axis results from PGN
        Second_dataset = Pandora;      
        Second_instrument_name = 'Pandora';
        binwidth = 2;
        Second_dataset.target2 = Second_dataset.NO2_surf;% the PGN multi-axis data product
        Env.target = Env.vmr(:,1);% EnviMes surface mixing ratio (i.e., first layer)
        measurand_txt = 'Surface NO_2 [ppbv]';% absolute observation
        measurand_txt2 = '(Surface NO_2 ratio)';% ratio text      
        measurand_txt_target2 = measurand_txt;
        avg_time = 10; % [min] EnviMes and Pandora data will be averged into this time bin        
        yrange = [0 40];
    case 'in_situ'
        load('C:\Projects\surface_NO2\from_DrSu\formated\MECP_Stn_34021.mat');% NAPS in situ observations at Donsview
        Second_dataset = MECP;      
        Second_instrument_name = 'In situ';
        binwidth = 2;
        Second_dataset.target2 = Second_dataset.NO2;% in situ NO2 results
        Env.target = Env.vmr(:,1);% EnviMes surface mixing ratio (i.e., first layer)
        measurand_txt = 'Surface NO_2 [ppbv]';% absolute observation
        measurand_txt2 = '(Surface NO_2 ratio)';% ratio text    
        measurand_txt_target2 = measurand_txt;
        avg_time = 60; % [min] EnviMes and in situ data will be averged into this time bin    
        yrange = [0 40];
    case 'in_situ_PM'% testing only!!
        load('C:\Projects\surface_NO2\from_DrSu\formated\MECP_Stn_34021.mat');% NAPS in situ observations at Donsview
        use_yyaxis = true;
        Second_dataset = MECP;      
        Second_instrument_name = 'In situ';
        binwidth = 'auto';
        Second_dataset.target2 = Second_dataset.PM;% in situ PM2.5/20 --> note this is just for testing! 
%         Env.target = Env.tvcd(:,1);% EnviMes AOD
        Env.target = Env.profile(:,1);% EnviMes AOD
        measurand_txt = 'Surface Extinction [km^-^1]';% absolute observation
        measurand_txt2 = '(SE/PM2.5 ratio)';% ratio text    
        measurand_txt_target2 = 'PM 2.5 [\mug/m^3]';
        avg_time = 10; % [min] EnviMes and in situ data will be averged into this time bin   
        yrange = [0 2];
    case 'AEROCAN'
        load('C:\Projects\AEROCAN\output\GTA_SDA_V3L15-2013_2021.mat');% this is the SDA version, re-processed by Ihab with different wavelength
        tf_site= strcmp(data.AERONET_Site,'Toronto');
        data(~tf_site,:) = [];
        data.AERONET_Site = [];
        data.Dateddmmyyyy = [];      
        data.Timehhmmss = [];   
        Second_dataset = data;
        Second_instrument_name = 'AEROCAN';
        binwidth = 0.1;
        Second_dataset.target2 = Second_dataset.Total_AOD_500nmtau_a;% SDA total AOD
        Env.target = Env.tvcd(:,1);% EnviMes AOD
        measurand_txt = 'AOD';% absolute observation
        measurand_txt2 = '(AOD ratio)';% ratio text    
        measurand_txt_target2 = measurand_txt;
        avg_time = 10; % [min] EnviMes and AERODCAN data will be averged into this time bin     
        
        tf = Second_dataset.target2 < 0;% simple QC
        Second_dataset(tf,:) = [];
        yrange = [0 1];
    case 'AEROCAN interp1'
        data_path = 'C:\Projects\AEROCAN\output\';
        data = [];
        for i =2016:2019
            load([data_path 'Toronto_AOD_V3L20-1999_2020_EnviMes_wv_' num2str(i) '.mat']);% the path for interped AOD
            data = [data;output];
        end
        tf_site= strcmp(data.AERONET_Site,'Toronto');
        data(~tf_site,:) = [];
        data.AERONET_Site = [];
        data.Dateddmmyyyy = [];      
        data.Timehhmmss = [];   
        data.AERONET_Site_Name = [];
        Second_dataset = data;
        Second_instrument_name = 'AEROCAN';
        binwidth = 0.1;
        Second_dataset.target2 = Second_dataset.AOD_target_1;% interped AOD 360 nm 
        Env.target = Env.tvcd(:,1);% EnviMes AOD
        measurand_txt = 'AOD';% absolute observation
        measurand_txt2 = '(AOD ratio)';% ratio text    
        measurand_txt_target2 = measurand_txt;
        avg_time = 10; % [min] EnviMes and AERODCAN data will be averged into this time bin     
        
        tf = Second_dataset.target2 < 0;% simple QC
        Second_dataset(tf,:) = [];
        yrange = [0 1];        
    case 'AEROCAN interp2'
        data_path = 'C:\Projects\AEROCAN\output\';
        data = [];
        for i =2016:2019
            load([data_path 'Toronto_AOD_V3L20-1999_2020_EnviMes_wv_' num2str(i) '.mat']);% the path for interped AOD
            data = [data;output];
        end
        tf_site= strcmp(data.AERONET_Site,'Toronto');
        data(~tf_site,:) = [];
        data.AERONET_Site = [];
        data.Dateddmmyyyy = [];      
        data.Timehhmmss = [];   
        data.AERONET_Site_Name = [];
        Second_dataset = data;
        Second_instrument_name = 'AEROCAN';
        binwidth = 0.1;
        Second_dataset.target2 = Second_dataset.AOD_target_2;% interped AOD 477 nm 
        Env.target = Env.tvcd(:,1);% EnviMes AOD
%         Env.target = Env.profile(:,1)
        measurand_txt = 'AOD';% absolute observation
        measurand_txt2 = '(AOD ratio)';% ratio text    
        measurand_txt_target2 = measurand_txt;
        avg_time = 10; % [min] EnviMes and AERODCAN data will be averged into this time bin     
        
        tf = Second_dataset.target2 < 0;% simple QC
        Second_dataset(tf,:) = [];
        yrange = [0 1];         
    case 'FTIR'
        data_path = 'D:\Projects\FTIR\reformat\';
        load([data_path 'FTIR_HCHO.mat']);% TAO FTIR data
        Second_dataset = data;
        Second_instrument_name = 'FTIR';
        binwidth = 0.1;
        Second_dataset.target2 = Second_dataset.int_4km_partial_column;% 0-4 km integrated column [DU]
        Env.target = Env.tvcd(:,1)./DU;% EnviMes HCHO column
        measurand_txt = 'HCHO [DU]';% absolute observation
        measurand_txt2 = '(HCHO ratio)';% ratio text    
        measurand_txt_target2 = measurand_txt;
        avg_time = 30; % [min] EnviMes and FTIR data will be averged into this time bin     
        yrange = [0 1];                 
    case 'FTIR_surface'
        data_path = 'D:\Projects\FTIR\reformat\';
        load([data_path 'FTIR_HCHO.mat']);% TAO FTIR data
        Second_dataset = data;
        Second_instrument_name = 'FTIR';
        binwidth = 0.1;
        Second_dataset.target2 = Second_dataset.vmr(:,end).*1e3;% convert vmr from ppmv to ppbv
        Env.target = Env.vmr(:,1);% EnviMes HCHO vmr in ppbv
        measurand_txt = 'HCHO [ppbv]';% absolute observation
        measurand_txt2 = '(HCHO ratio)';% ratio text    
        measurand_txt_target2 = measurand_txt;
        avg_time = 30; % [min] EnviMes and FTIR data will be averged into this time bin     
        yrange = [0 10];           
end

%% check sampling interval, if needed
Env_dummy.UTC = Env.UTC;
Env_dummy.UTC(1,:) = [];
Env_dummy.UTC(end+1,:) = Env_dummy.UTC(end,:) + minutes(5);
dt = Env_dummy.UTC - Env.UTC;
tf = (dt > minutes(30)) | (dt < -minutes(30));
dt(tf,:) = [];
figure;
subplot(1,2,1);hold all;
histogram(dt);
ylabel([instrument_name ' f']);
xlabel([instrument_name ' sampling interval']);

P_dummy.UTC = Second_dataset.UTC;
P_dummy.UTC(1,:) = [];
%P_dummy.UTC(end+1,:) = P_dummy.UTC(end,:) + minutes(5);
P_dummy.UTC(end+1,:) = NaT;
dt = P_dummy.UTC - Second_dataset.UTC;
if strcmp(comparison_with,'in_situ')
    disp('Comparison with in situ observations');
    tf = (dt > minutes(60)) | (dt < -minutes(60));
    dt(tf,:) = [];
else
    disp('Comparison with high-temporal resolution data, remove observations interval > 30 min in histogram');
    tf = (dt > minutes(30)) | (dt < -minutes(30));
    dt(tf,:) = [];
end

subplot(1,2,2);hold all;
histogram(dt);
ylabel([Second_instrument_name ' f']);
xlabel([Second_instrument_name ' sampling interval']);

fig_name = [instrument_name '_' target '_vs_' Second_instrument_name '_' comparison_with '_sampling_interval'];
print_setting(1/2,save_fig,[plot_path fig_name]);



%% merging data
min_year = min([min(Second_dataset.UTC.Year) min(Env.UTC.Year)]);% get the starting year
max_year = max([max(Second_dataset.UTC.Year) max(Env.UTC.Year)]);% get the ending year
newTimes  = [datetime(min_year,1,1):minutes(avg_time):datetime(max_year,12,31)];
if ~istimetable(Second_dataset)
    Second_dataset = table2timetable(Second_dataset);
end
Env = table2timetable(Env);
Second_dataset_retime = retime(Second_dataset,newTimes ,'mean');
Env_retime = retime(Env,newTimes,'mean');

C = innerjoin(Env_retime,Second_dataset_retime,'RightKeys','UTC','LeftKeys','UTC');

%% simple QAQC
tf = isnan(C.tvcd);
C(tf,:) = [];

tfa = C.UTC.Hour-5 <=0;
C(tfa,:) = [];

tf = isnan(C.target) | isnan(C.target2);
C(tf,:) = [];

%% time serise
figure;
subplot(1,2,1);hold all;
if use_yyaxis
    yyaxis left
    plot(Env.UTC,Env.target,'.');% EnviMes
    ylabel(measurand_txt);
    yyaxis right
    plot(Second_dataset.UTC,Second_dataset.target2,'.');% e.g., Pandora 
    ylabel(measurand_txt_target2);
else
    plot(Env.UTC,Env.target,'.');% EnviMes
    plot(Second_dataset.UTC,Second_dataset.target2,'.');% e.g., Pandora 
    ylabel(measurand_txt);
end
xlabel('Date');
legend({instrument_name,Second_instrument_name});
title('Raw data');
subplot(1,2,2);hold all;
if use_yyaxis
    yyaxis left
    plot(C.UTC,C.target,'.');% EnviMes
    ylabel(measurand_txt);
    yyaxis right
    plot(C.UTC,C.target2,'.');% e.g., Pandora 
    ylabel(measurand_txt_target2);
else
    plot(C.UTC,C.target,'.');% EnviMes
    plot(C.UTC,C.target2,'.');% e.g., Pandora 
    ylabel(measurand_txt);
end
xlabel('Date');
legend({instrument_name,Second_instrument_name});
title('Coincident data');
fig_name = [instrument_name '_' target '_vs_' Second_instrument_name '_' comparison_with '_timeserise'];
print_setting(1/2,save_fig,[plot_path fig_name]);

%% comparison


% simple scatter for trop column
figure;
if use_yyaxis
    subplot(2,2,[1 3]);hold all;
    line_fits_local(C.target2,C.target,false,false);% e.g., Pandora vs. EnviMes
%     x_min = min([min(C.target2) min(C.target) 0]);
%     x_max = max([max(C.target2) max(C.target)]);
%     xlim([x_min x_max]);
%     ylim([x_min x_max]);
    xlabel([Second_instrument_name ' ' measurand_txt_target2]);
    ylabel([instrument_name ' ' measurand_txt]);
    subplot(2,2,2);hold all;
    histogram(C.target,'Normalization','probability');
    legend({instrument_name});
    ylabel('f');
    xlabel(measurand_txt);
    subplot(2,2,4);hold all;
    histogram(C.target2,'Normalization','probability'); 
    legend({Second_instrument_name});
    ylabel('f');
    xlabel(measurand_txt_target2);    
else
    subplot(1,2,1);hold all;
    line_fits_local(C.target2,C.target);% e.g., Pandora vs. EnviMes
    x_min = min([min(C.target2) min(C.target) 0]);
    x_max = max([max(C.target2) max(C.target)]);
    xlim([x_min x_max]);
    ylim([x_min x_max]);
    xlabel([Second_instrument_name ' ' measurand_txt_target2]);
    ylabel([instrument_name ' ' measurand_txt]);
    subplot(1,2,2);hold all;
    histogram(C.target,'Normalization','probability','BinWidth',binwidth);
    histogram(C.target2,'Normalization','probability','BinWidth',binwidth);
    legend({instrument_name,Second_instrument_name});
    ylabel('f');
    xlabel(measurand_txt);
end

fig_name = [instrument_name '_' target '_vs_' Second_instrument_name '_' comparison_with '_scatter_and_hist'];
print_setting(1/2,save_fig,[plot_path fig_name]);

% seasonal difference for trop column with boxplot
tf = abs(C.target./C.target2) < 10;% a simple filter for ratio plot
figure;
subplot(1,2,1);hold all;
plot(datetime(0,C.UTC.Month(tf,:),C.UTC.Day(tf,:)), C.target(tf,:)./C.target2(tf,:),'.');% e.g., EnviMes/Pandora
xlabel(['Date']);
ylabel([instrument_name '/' Second_instrument_name ' ' measurand_txt2]);
subplot(1,2,2);hold all;
boxplot(C.target(tf,:)./C.target2(tf,:),C.UTC.Month(tf,:)) ;
xlabel(['Month']);
ylabel([instrument_name '/' Second_instrument_name ' ' measurand_txt2]);
fig_name = [instrument_name '_' target '_vs_' Second_instrument_name '_' comparison_with '_seasonal'];
print_setting(1/2,save_fig,[plot_path fig_name]);



% diurnal difference for trop column with boxplot
figure;
subplot(1,2,1);hold all;
plot(C.sza(tf,:), C.target(tf,:)./C.target2(tf,:),'.');% e.g., EnviMes/Pandora
xlabel(['SZA [degrees]']);
ylabel([instrument_name '/' Second_instrument_name  ' ' measurand_txt2]);
subplot(1,2,2);hold all;
boxplot(  C.target(tf,:)./C.target2(tf,:),round(C.sza(tf,:).*10,-2)./10);% EnviMes/Pandora
xlabel(['SZA [degrees]']);
ylabel([instrument_name '/' Second_instrument_name  ' ' measurand_txt2]);
fig_name = [instrument_name '_' target '_vs_' Second_instrument_name '_' comparison_with '_diurnal'];
print_setting(1/2,save_fig,[plot_path fig_name]);

% seasonal and diurnal with boxplot
if use_yyaxis
    c_map = lines(2);
    figure;
    subplot(1,3,1);hold all;% seasonal 
    X1 = 1.3:12.3;
    yyaxis left;
    boxplot(C.target(tf,:),C.UTC.Month(tf,:),'positions', X1, 'labels', X1,'Widths',0.2,'Colors',c_map(1,:),'Symbol','') ;
    m1 = grpstats(C.target(tf,:),C.UTC.Month(tf,:),'median');
    plot(X1,m1,'-','Color',c_map(1,:));
    ylabel([measurand_txt]);
    X2 = 1:12;
    yyaxis right;
    boxplot(C.target2(tf,:),C.UTC.Month(tf,:),'positions', X2, 'labels', X2,'Widths',0.2,'Colors',c_map(2,:),'Symbol','') ;
    m2 = grpstats(C.target2(tf,:),C.UTC.Month(tf,:),'median');
    plot(X2,m2,'-','Color',c_map(2,:));
    xlabel(['Month']);
    ylabel([measurand_txt_target2]);
    legend({[instrument_name ' [' strrep(target,'_',' ') ']'],[Second_instrument_name ' [' strrep(comparison_with,'_',' ') ']']});
%     ylim(yrange);

    subplot(1,3,2);hold all;% diurnal
    X1 = 22:10:82;
    yyaxis left;
    boxplot(C.target(tf,:),round(C.sza(tf,:).*10,-2)./10,'positions', X1, 'labels', X1,'Widths',1,'Colors',c_map(1,:),'Symbol','') ;    
    m1 = grpstats(C.target(tf,:),round(C.sza(tf,:).*10,-2)./10,'median');
    plot(X1,m1,'-','Color',c_map(1,:));
    ylabel([measurand_txt]);
    X2 = 20:10:80;
    yyaxis right;
    boxplot(C.target2(tf,:),round(C.sza(tf,:).*10,-2)./10,'positions', X2, 'labels', X2,'Widths',1,'Colors',c_map(2,:),'Symbol','') ;
    xlabel(['SZA [degrees]']);    
    m2 = grpstats(C.target2(tf,:),round(C.sza(tf,:).*10,-2)./10,'median');    
    plot(X2,m2,'-','Color',c_map(2,:));
    ylabel([measurand_txt_target2]);
    legend({[instrument_name ' [' strrep(target,'_',' ') ']'],[Second_instrument_name ' [' strrep(comparison_with,'_',' ') ']']});
%     ylim(yrange);

    subplot(1,3,3);hold all;% diurnal 2
    obs_hours = unique(C.UTC.Hour(tf,:)-5);
    % X1 = min(obs_hours)+0.3:max(obs_hours)+0.3;
    X1 = double(obs_hours)+0.3;
    yyaxis left;
    boxplot(C.target(tf,:),C.UTC.Hour(tf,:)-5,'positions', X1, 'labels', X1,'Widths',0.2,'Colors',c_map(1,:),'Symbol','') ;
    m1 = grpstats(C.target(tf,:),C.UTC.Hour(tf,:)-5,'median');
    plot(X1,m1,'-','Color',c_map(1,:));
    ylabel([measurand_txt]);
    X2 = double(obs_hours);
    yyaxis right;
    boxplot(C.target2(tf,:),C.UTC.Hour(tf,:)-5,'positions', X2, 'labels', X2,'Widths',0.2,'Colors',c_map(2,:),'Symbol','') ;    
    m2 = grpstats(C.target2(tf,:),C.UTC.Hour(tf,:)-5,'median');    
    plot(X2,m2,'-','Color',c_map(2,:));
    xlabel(['LST [hour]']);
    ylabel([measurand_txt_target2]);
    legend({[instrument_name ' [' strrep(target,'_',' ') ']'],[Second_instrument_name ' [' strrep(comparison_with,'_',' ') ']']});
%     ylim(yrange);
else
    c_map = lines(2);
    figure;
    subplot(1,3,1);hold all;% seasonal 
    X1 = 1.3:12.3;
    boxplot(C.target(tf,:),C.UTC.Month(tf,:),'positions', X1, 'labels', X1,'Widths',0.2,'Colors',c_map(1,:),'Symbol','') ;
    m1 = grpstats(C.target(tf,:),C.UTC.Month(tf,:),'median');
    plot(X1,m1,'-','Color',c_map(1,:));
    X2 = 1:12;
    boxplot(C.target2(tf,:),C.UTC.Month(tf,:),'positions', X2, 'labels', X2,'Widths',0.2,'Colors',c_map(2,:),'Symbol','') ;
    m2 = grpstats(C.target2(tf,:),C.UTC.Month(tf,:),'median');
    plot(X2,m2,'-','Color',c_map(2,:));
    xlabel(['Month']);
    ylabel([measurand_txt]);
    legend({[instrument_name ' [' strrep(target,'_',' ') ']'],[Second_instrument_name ' [' strrep(comparison_with,'_',' ') ']']});
    ylim(yrange);

    subplot(1,3,2);hold all;% diurnal
    X1 = 22:10:82;
    boxplot(C.target(tf,:),round(C.sza(tf,:).*10,-2)./10,'positions', X1, 'labels', X1,'Widths',1,'Colors',c_map(1,:),'Symbol','') ;
    X2 = 20:10:80;
    boxplot(C.target2(tf,:),round(C.sza(tf,:).*10,-2)./10,'positions', X2, 'labels', X2,'Widths',1,'Colors',c_map(2,:),'Symbol','') ;
    xlabel(['SZA [degrees]']);
    m1 = grpstats(C.target(tf,:),round(C.sza(tf,:).*10,-2)./10,'median');
    m2 = grpstats(C.target2(tf,:),round(C.sza(tf,:).*10,-2)./10,'median');
    plot(X1,m1,'-','Color',c_map(1,:));
    plot(X2,m2,'-','Color',c_map(2,:));
    ylabel([measurand_txt]);
    legend({[instrument_name ' [' strrep(target,'_',' ') ']'],[Second_instrument_name ' [' strrep(comparison_with,'_',' ') ']']});
    ylim(yrange);

    subplot(1,3,3);hold all;% diurnal 2
    obs_hours = unique(C.UTC.Hour(tf,:)-5);
    % X1 = min(obs_hours)+0.3:max(obs_hours)+0.3;
    X1 = double(obs_hours)+0.3;
    boxplot(C.target(tf,:),C.UTC.Hour(tf,:)-5,'positions', X1, 'labels', X1,'Widths',0.2,'Colors',c_map(1,:),'Symbol','') ;
    % X2 = min(obs_hours):max(obs_hours);
    X2 = double(obs_hours);
    boxplot(C.target2(tf,:),C.UTC.Hour(tf,:)-5,'positions', X2, 'labels', X2,'Widths',0.2,'Colors',c_map(2,:),'Symbol','') ;
    xlabel(['LST [hour]']);
    m1 = grpstats(C.target(tf,:),C.UTC.Hour(tf,:)-5,'median');
    m2 = grpstats(C.target2(tf,:),C.UTC.Hour(tf,:)-5,'median');
    plot(X1,m1,'-','Color',c_map(1,:));
    plot(X2,m2,'-','Color',c_map(2,:));
    ylabel([measurand_txt]);
    legend({[instrument_name ' [' strrep(target,'_',' ') ']'],[Second_instrument_name ' [' strrep(comparison_with,'_',' ') ']']});
    ylim(yrange);
end

fig_name = [instrument_name '_' target '_vs_' Second_instrument_name '_' comparison_with '_ssl'];
print_setting(1/3,save_fig,[plot_path fig_name]);


function [intercept,slope,slope_nlm,mdl_lm,mdl_nlm] = line_fits_local(x,y,only_simple,with_1_on_1)
% simple linear fit
if nargin ==2
    only_simple = false;
    with_1_on_1 = true;
end

mdl_lm = fitlm(x,y);

myfun = @(b,x)b(1)*x;
beta =[1];
mdl_nlm = fitnlm(x,y,myfun,beta);
% mdl_nlm = fitlm(x,y,'Intercept',false);% this is another way to fit with intercept set to 0
% figure;hold all;
plot(x,y,'.');
plot(x,predict(mdl_lm,x),'color',[0.1 0.1 0.9]);
if ~only_simple
    plot(x,predict(mdl_nlm,x),'color',[0.9 0.1 0.5]);
end
if with_1_on_1
    plot(x,x,'k--');
end
intercept = mdl_lm.Coefficients.Estimate(1);
slope = mdl_lm.Coefficients.Estimate(2);
slope_nlm = mdl_nlm.Coefficients.Estimate(1);
x0 = double((max(x) - min(x))/8);
y0 = double((max(y) - min(y))/5);
text(x0,y0*4, ['y = ' num2str(slope) '*x + ' num2str(intercept)],'color',[0.1 0.1 0.9]);
if ~only_simple
    text(x0,y0*3.5,['y = ' num2str(slope_nlm) '*x'],'color',[0.9 0.1 0.5]);
end
text(x0,y0*3,['N = ' num2str(numel(x)) ]);
text(x0,y0*2.5,['R = ' num2str(corr(x,y)) ]);

