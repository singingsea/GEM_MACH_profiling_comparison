function read_gem_mach_no2_profile()

% site = 'Downsview';
site = 'Toronto-West';
traget_gas = 'NO2';
read_files = true;% if true, the code will read data from "raw" nc files
data_folder = 'D:\Projects\GEM_MACH\profile_data\Toronto\';
output_folder = ['D:\Projects\GEM_MACH\profile_data_output\' site '\' ];mkdir(output_folder);
save_fig = 1;
plot_path = ['D:\Projects\GEM_MACH\profile_data_plots\' site '\' ];mkdir(plot_path);


switch site
    case 'Downsview'
        user_lat=43.7810; % Downsview
        user_lon=-79.4680;
    case 'Egbert'
        user_lat=44.2300; %  Egbert
        user_lon=-79.7800;
    case 'UTSG'
        user_lat=43.6605;% St George
        user_lon=-79.39860;
    case 'UTSC'        
        user_lat=43+47/60; % UTSC
        user_lon=-(79 + 11/60);
    case 'Toronto-West'        
        user_lat=43.7109; % Toronto West
        user_lon=-79.5418;
end

if read_files
    list = dir([data_folder '*.nc']);
    GEM_profile = table;
    for i_file = 1:numel(list)
        f = [data_folder list(i_file).name];% nc fie that will be read
        disp(['Reading file: ' f ]);
        timestamp = list(i_file).name(end-10:end-3);% yyyyMMdddouble check this if the file name is changed
        yyyy = str2num(timestamp(1:4));MM = str2num(timestamp(5:6));dd = str2num(timestamp(7:8));
        t_hours = ncread(f,'/hours');% 
        GZt= ncread(f,'/GZt');% The GZt is in decameters ASL
        PX= ncread(f,'/PX');% hPa
        TT= ncread(f,'/TT');% C
        TNO2= ncread(f,'/TNO2');% ppb
        lat= ncread(f,'/lat');% deg
        lon= ncread(f,'/lon');% deg
        VCD= ncread(f,'/VCD');% molec/cm2; 'Partial column up to 4km'
        UTC = datetime(yyyy,MM,dd,t_hours,0,0);% generate the UTC timestamp

        N = size(lat);
        for i = 1:N(1)% loop over lat grids
            for j = 1:N(2)% loop over lon grids
                lat_1grid = lat(i,j);
                lon_1grid = lon(i,j);        
                d(i,j) = get_distance(user_lat,user_lon,lat_1grid,lon_1grid);% calculate the distance between model grid to site
            end
        end
        closest_d = min(d(:));
        [lat_idx,lon_idx] = find(d == closest_d);
        min_lat = lat(lat_idx,lon_idx);
        min_lon = lon(lat_idx,lon_idx)-360;
        disp(['Site: ' site ' @ lat:' num2str(user_lat) ', lon:'  num2str(user_lon)]);
        disp(['Closest GEM-MACH grid centre: ' num2str(closest_d) 'km @ lat:' num2str(min_lat) ', lon:'  num2str(min_lon)]);

        % extract the point that closest to site
        h = squeeze(GZt(lat_idx,lon_idx,:,:))./100;% convert to height in km [lat, lon, level, time]
        P = squeeze(PX(lat_idx,lon_idx,:,:));% in hPa
        T = squeeze(TT(lat_idx,lon_idx,:,:));% in C
        NO2_vmr = squeeze(TNO2(lat_idx,lon_idx,:,:));% in ppb
        NO2_VCD = squeeze(VCD(lat_idx,lon_idx,:,:));% in molec/cm2

        GEM_profile.UTC = UTC;
        GEM_profile.h = h';
        GEM_profile.VCD = NO2_VCD;
        GEM_profile.NO2_vmr = NO2_vmr';
        GEM_profile.P = P';
        GEM_profile.T = T';

        if i_file == 1
            GEM_profiles = GEM_profile;
        else
            GEM_profiles = [GEM_profiles;GEM_profile];
        end

    end

    save([output_folder 'GEM_NO2_profiles_' site],'GEM_profiles');
end

%% plotting
if ~read_files
    load(['D:\Projects\GEM_MACH\profile_data_output\Downsview\GEM_NO2_profiles_' site '.mat']);
end
data = GEM_profiles;
tf = data.NO2_vmr > 500;% simple QC
tf_bad_data = logical(sum(tf,2));
data(tf_bad_data,:) = [];


interp_timestamp = min(data.UTC):minutes(60):max(data.UTC);
interp_method = 'mean';
% interp_method = 'linear';
data = table2timetable(data);
data = retime(data,interp_timestamp,interp_method);
tf_nan = sum(isnan(data.NO2_vmr),2) > 0;
data(tf_nan,:) = [];
data.UTC_datenum = datenum(data.UTC);
h_interp = 0:0.02:4;
X = data.UTC_datenum';
Y = data.h;
V = data.NO2_vmr';
Xq= X;
Yq = h_interp;
% Vq = interp1(Y(1,:),V,Yq);
for i = 1:height(data)
        try
        Vq(:,i) = interp1(Y(i,:),V(:,i),Yq);
        catch 
            a
        end
end
data.LST = data.UTC - hours(5);



for i = min(data.UTC.Year):max(data.UTC.Year)% loop over years
    figure;
    for j = 1:12 % loop over months
        tf_1yr_1month = (data.UTC.Year == i) & (data.UTC.Month == j);
        data_1yr_1month = data(tf_1yr_1month,:);
        subplot(12,1,j);hold all;
        if j == 1
            title(['GEM-MACH ' traget_gas ' profiles @ ' site ' ' num2str(i)]);
        end
        xlim([datenum(datetime(i,j,1)) datenum(datetime(i,j+1,1))]);datetick('x','mmm-dd','keeplimits')
        %         datetick('x','ddd','keeplimits')
%         ylim([0 4]);
        ylim([0 1.5]);
        caxis([0 45]);
        if ~isempty(data_1yr_1month)
            plot_profile(data.UTC(tf_1yr_1month,:),h_interp,Vq(:,tf_1yr_1month), ' ');
            if j == 6
                ylabel('Altitude [km]');
            end
            
            if j == 6
                hcb = colorbar('Location','manual');
                hcb.Position = [0.93, 0.1, 0.02, 0.8];%  [left, bottom, width, height]
                colorTitleHandle = get(hcb,'Title');
                titleString = [traget_gas ' [ppbv]'];
                current_position = colorTitleHandle.Position;
                current_position(1,1) = -10;
                current_position(1,2) = 550;
                set(colorTitleHandle ,'String',titleString,'Position',current_position);
            end
        end
        
    end
    fig_name = ['GEM_MACH_' traget_gas '_profiles_' num2str(i)];
    print_setting(1,save_fig,[plot_path fig_name]);
    
end
% figure;hold all;
% plot(NO2_vmr(:,:),h(:,:));



function plot_profile(UTC,h,profile,colorTitle)
mycolormap = customcolormap(linspace(0,1,12), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691','#ffffff',});
% mycolormap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'});
colormap(mycolormap);


imagesc(datenum(UTC),h,profile);
% contourf(datenum(UTC),h,profile) 
set(gca, 'YDir','normal');
% datetick('x','mmm-dd','keeplimits')
% datetick('x','yyyy','keeplimits')
% hcb = colorbar;
% colorTitleHandle = get(hcb,'Title');
% titleString = colorTitle;
% set(colorTitleHandle ,'String',titleString);
% ylabel('Altitude [m]');
grid on;