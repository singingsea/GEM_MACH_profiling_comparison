function read_NDACC_MAX()
% this function is to read NDACC centralized processed MAX-DOAS data
addpath(genpath('C:\Users\ZhaoX\Documents\MATLAB\matlab'));

% instrument = '1694_1';
instrument = '1694_2';
years = 2016:2021;
% years = 2021;

for i = 1:numel(years)   
    year = num2str(years(i));
    switch instrument
        case '1694_2'
            data_folder = ['D:\Projects\EnviMes\NDACC_processed\tropo\stations\TORONTO.DOWNSVIEW\' instrument '\o4-no2,o4-hcho\fv002\' year '\'];
            tragets = {'NO2_428NM','O4_360NM','HCHO_343NM'};
            
        case '1694_1'
            data_folder = ['D:\Projects\EnviMes\NDACC_processed\tropo\stations\TORONTO.DOWNSVIEW\' instrument '\o4-no2\fv002\' year '\'];
            tragets = {'NO2_460NM','O4_477NM'};
    end

    output_path = ['D:\Projects\EnviMes\output_profiles\' instrument '\'];
%     plot_path = [output_path 'plot\' instrument '\'];                
    mkdir(output_path);
%     mkdir(plot_path);
    
    for j = 1:numel(tragets)
        target = char(tragets(j));
        disp(string(['Reading: EnviMes' instrument ' ' target ' in ' year]));
        read_one_target(instrument,year,data_folder,target,output_path);
    end
end
% target = 'NO2_428NM';% 
% target = 'O4_360NM';% 
% target = 'HCHO_343NM';% 

function read_one_target(instrument,year,data_folder,target,output_path)

list = dir([data_folder '*.nc']);
files = list(~[list.isdir],:);

UTC_all = [];
profile_all = [];
profile_error_all = [];
tvcd_all = [];
tvcd_error_all = [];
qa_all = [];
dof_all = [];
sza_all = [];
saz_all = [];
% h_all = [];
thickness_of_layer_all = [];
air_in_column_all = [];


for i = 1:numel(files)    
    file = [data_folder files(i).name];
    disp(['Reading file: ' files(i).name]);
    try
        if contains(target,'O4')
            profile = ncread(file,['/PROFILE/MMF/' target '/aerosol_extinction_profile']);% 'retrieved aerosol profile' ['km-1']
            profile_error = ncread(file,['/PROFILE/MMF/' target '/aerosol_extinction_profile_error']);% 'combined smoothing and noise error of retrieved aerosol profile' ['km-1']
            tvcd = ncread(file,['/PROFILE/MMF/' target '/aerosol_optical_depth']);% 'Total aerosol optical depth' 
            tvcd_error = ncread(file,['/PROFILE/MMF/' target '/aerosol_optical_depth_error']);% 'Total aerosol optical depth error calculated from covariance measurement noise matrix and covariance smoothing error matrix' 
        else
            profile = ncread(file,['/PROFILE/MMF/' target '/concentration_profile']);% read the retreived profile ['molec/cm3']
            profile_error = ncread(file,['/PROFILE/MMF/' target '/concentration_profile_error']);% combined smoothing and noise error of retrieved gas profile ['molec/cm3']
            tvcd = ncread(file,['/PROFILE/MMF/' target '/tropospheric_vertical_column_density']);% read the trop column ['molec/cm2']
            tvcd_error = ncread(file,['/PROFILE/MMF/' target '/tropospheric_vertical_column_density_error']);% read the trop column error ['molec/cm2']
        end
        
        qa = ncread(file,['/PROFILE/MMF/' target '/qa_flag']);% read qa values
        dof = ncread(file,['/PROFILE/MMF/' target '/degrees_of_freedom']);% degree of freedom    

        sza = ncread(file,'/PROFILE/solar_zenith_angle');
        saz = ncread(file,'/PROFILE/solar_azimuth_angle');        
        yyyy = ncread(file,'/PROFILE/year');
        fday = ncread(file,'/PROFILE/fractional_day_of_year');
        h = ncread(file,'/PROFILE/height_of_layer_center_above_sea_level');
        thickness_of_layer = ncread(file,'/PROFILE/thickness_of_layer');%'thickness of retrieval layer' [m]
        n_dim = size(sza);thickness_of_layer = repmat(thickness_of_layer,1,n_dim(1));% expand thinknees of layers to same dim of other parameters
        UTC = datetime(yyyy,1,1) + days(fday - 1);    
        air_in_column = ncread(file,'/PROFILE/MMF/air_column_in_layer');%'partial column density of air in layer' ['molec/cm2']
        
        UTC_all = [UTC_all;UTC];
        profile_all = [profile_all,profile];
        profile_error_all = [profile_error_all,profile_error];
        tvcd_all = [tvcd_all;tvcd];
        tvcd_error_all = [tvcd_error_all;tvcd_error];
        qa_all = [qa_all;qa];
        dof_all = [dof_all;dof];
        sza_all = [sza_all;sza];
        saz_all = [saz_all;saz];
    %     h_all = [h_all;h];
        thickness_of_layer_all = [thickness_of_layer_all,thickness_of_layer];
        air_in_column_all = [air_in_column_all,air_in_column];
    
    
    catch
        disp(['Failed to read file: ' file]);
    end

end
data = table;
data.UTC = UTC_all;
data.tvcd = tvcd_all;
data.tvcd_error = tvcd_error_all;
data.qa = qa_all;
data.dof = dof_all;
data.sza = sza_all;
data.saz = saz_all;
data.profile = profile_all';
data.profile_error = profile_error_all';
data.thickness_of_layer = thickness_of_layer_all';
data.air_in_column = air_in_column_all';
data.vmr = data.profile./(data.air_in_column./(data.thickness_of_layer.*100)).*1e9;% calculate vmr profile and convert to ppbv

% save data
save([output_path target '_' instrument '_' year '.mat'],'data','h','instrument','year','data_folder','target');

% 
% interp_timestamp = 'hourly';
% interp_method = 'linear';
% data = table2timetable(data);
% data = retime(data,interp_timestamp,interp_method);

% data.UTC_datenum = datenum(data.UTC);
% h_interp = 0:1:4000;
% X = data.UTC_datenum';
% Y = h';
% V = data.profile';
% Xq= X;
% Yq = h_interp;
% Vq = interp1(Y,V,Yq);
% 
% 
% figure;hold all;
% plot_profile(data.UTC,h,data.profile',['NO_2 [molec/cm^3]']);
% ylim([0 4000]);
% fig_name = [target '_' instrument '_' year '_raw'];
% print_setting(1/2,save_fig,[plot_path fig_name]);
% 
% 
% figure;hold all;
% plot_profile(data.UTC,h_interp,Vq,['NO_2 [molec/cm^3]']);
% ylim([300 4000]);
% x_min = min(datenum(data.UTC));x_max = max(datenum(data.UTC));
% xlim([x_min x_max]);
% fig_name = [target '_' instrument '_' year '_interp'];
% print_setting(1/2,save_fig,[plot_path fig_name]);

% %% plot
% function plot_profile(UTC,h,profile,colorTitle)
% mycolormap = customcolormap(linspace(0,1,12), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691','#ffffff',});
% % mycolormap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'});
% colormap(mycolormap);
% 
% 
% imagesc(datenum(UTC),h,profile);
% % contourf(datenum(UTC),h,profile) 
% set(gca, 'YDir','normal');
% datetick('x','mmm-dd','keeplimits')
% hcb = colorbar;
% colorTitleHandle = get(hcb,'Title');
% titleString = colorTitle;
% set(colorTitleHandle ,'String',titleString);
% ylabel('Altitude [m]');
% grid on;