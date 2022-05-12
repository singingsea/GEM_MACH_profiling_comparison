function EnviMes_vs_GEM_MACH()
% this function is a simple one to check the NO2 profile results between
% EnviMes and GEM-MACH.
% the GEM-MACH reusluts are from code "D:\Projects\GEM_MACH\code\read_gem_mach_no2_profile.m"
plot_path = 'D:\Projects\EnviMes\plot\GEM-MACH\';
save_fig = 1;

%% load GEM-MACH NO2 profiles
site = 'Downsview';load(['D:\Projects\GEM_MACH\profile_data_output\Downsview\GEM_NO2_profiles_' site '.mat']);

tf_nan = (sum(GEM_profiles.NO2_vmr > 1e3,2)) > 0;
GEM_profiles(tf_nan,:) = [];

GEM_profiles.h = GEM_profiles.h - 0.187; % Location altitude [m]: 187; convert ASL to AGL by subtract Pandora height?

%% load EnviMes profiles

target = 'NO2_428NM';instrument = '1694_2';
% target = 'NO2_460NM';instrument = '1694_1';

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

%% merge profiles
data1 = table2timetable(GEM_profiles);% hourly
data2 = table2timetable(Env);% about 8 min 
data2 = retime(data2,'hourly','mean');% retime MAX-DOAS to hourly mean
tf_nan = isnan(data2.tvcd);
data2(tf_nan,:) = [];
data = innerjoin(data1,data2);

%% compare partial columns

tf_warm = (data.UTC.Month>=5) & (data.UTC.Month<=10);

figure;
fig_name = ['VCD_and_surface_comparison_' target];
print_setting(1,0,[plot_path fig_name]);

x = data.tvcd;% MAX-DOAS partial column
y = data.VCD;% GEM-MACH partial column
subplot(2,2,1);hold all;grid on;
line_fits_local2(x(tf_warm,:),y(tf_warm,:));
xlabel('MAX-DOAS NO_2 pVCD [molec/cm^2]');
ylabel('GEM-MACH NO_2 pVCD [molec/cm^2]');
title('Warm seasons');
xlim([0 1e17]);ylim([0 1e17]);
subplot(2,2,2);hold all;grid on;
line_fits_local2(x(~tf_warm,:),y(~tf_warm,:));
xlabel('MAX-DOAS NO_2 pVCD [molec/cm^2]');
ylabel('GEM-MACH NO_2 pVCD [molec/cm^2]');
title('Cold seasons');
xlim([0 1e17]);ylim([0 1e17]);

x = data.vmr(:,1);% MAX-DOAS partial column
y = data.NO2_vmr(:,end);% GEM-MACH bottom layer
subplot(2,2,3);hold all;grid on;
line_fits_local2(x(tf_warm,:),y(tf_warm,:));
xlabel('MAX-DOAS bottom layer VMR [ppbv]');
ylabel('GEM-MACH bottom layer VMR [ppbv]');
title('Warm seasons');
xlim([0 120]);ylim([0 120]);

subplot(2,2,4);hold all;grid on;
line_fits_local2(x(~tf_warm,:),y(~tf_warm,:));
xlabel('MAX-DOAS bottom layer VMR [ppbv]');
ylabel('GEM-MACH bottom layer VMR [ppbv]');
title('Cold seasons');
xlim([0 120]);ylim([0 120]);


print_setting(1,save_fig,[plot_path fig_name]);

%% regrid profiles to same vertical grids
h_interp = 0:0.2:4;% new grid in [km]

Y = data.h; % GEM-MACH profile regrid
V = data.NO2_vmr';% GEM-MACH NO2 profile
Yq = h_interp;
for i = 1:height(data)
        Vq(i,:) = interp1(Y(i,:),V(:,i),Yq);
end
Vq = Vq';

Y2 = 0.1:0.2:4; % MAX-DOAS profile regrid, 200 m resolution, 0-4 km
V2 = data.vmr';% MAX-DOAS NO2 profile
Yq = h_interp;
Vq2 = interp1(Y2,V2,Yq);

%% compare

%% profile plots
figure;
mycolormap = customcolormap(linspace(0,1,12), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691','#ffffff',});
subplot(2,1,1);
plot_profile(data.UTC,h_interp,Vq,mycolormap);
% ylim([0 1.5]);
colorbar;title('GEM-MACH NO_2 profile [ppbv]');
ylabel('Altitude [km]');
datetick('x','yyyy-mm','keeplimits');
subplot(2,1,2);
plot_profile(data.UTC,h_interp,Vq2,mycolormap);
% ylim([0 1.5]);
colorbar;title('MAX-DOAS NO_2 profile [ppbv]');
ylabel('Altitude [km]');
datetick('x','yyyy-mm','keeplimits');

fig_name = ['Profile_timeserise_' target];
print_setting(1/2,save_fig,[plot_path fig_name]);

%% calculate difference 
tf_warm = (data.UTC.Month>=5) & (data.UTC.Month<=10);% get cold and warm seasons index
Vq_warm = Vq(:,tf_warm);% vmr profiles in warm seasons
Vq_cold = Vq(:,~tf_warm);% vmr profiles in cold seasons
dVq = Vq - Vq2;% calculate difference
dVq_warm = dVq(:,tf_warm);% difference in warm seasons
dVq_cold = dVq(:,~tf_warm);% difference in cold seasons

dVq_mean = mean(dVq,2,'omitnan')';% calculate mean 
dVq_std = std(dVq,0,2,'omitnan')';% calcualte std
dVq_median = median(dVq,2,'omitnan')';% calculate median
tf_nan= isnan(dVq_mean);
p_dVq = (Vq - Vq2)./((Vq+Vq2)./2).*100;% % difference
% p_dVq = (Vq - Vq2)./Vq2.*100;
p_dVq_mean = mean(p_dVq,2,'omitnan')';
p_dVq_std = std(p_dVq,0,2,'omitnan')';
p_dVq_median = median(p_dVq,2,'omitnan')';

dVq_warm_mean = mean(dVq_warm,2,'omitnan')';% calculate mean in warm sesons
dVq_warm_std = std(dVq_warm,0,2,'omitnan')';% calcualte std
dVq_warm_median = median(dVq_warm,2,'omitnan')';% calculate median
tf_nan= isnan(dVq_warm_mean);
p_dVq_warm = (Vq - Vq2)./((Vq+Vq2)./2).*100;% % difference
% p_dVq_warm = (Vq - Vq2)./Vq2.*100;
p_dVq_warm_mean = mean(p_dVq_warm,2,'omitnan')';
p_dVq_warm_std = std(p_dVq_warm,0,2,'omitnan')';
p_dVq_warm_median = median(p_dVq_warm,2,'omitnan')';

dVq_cold_mean = mean(dVq_cold,2,'omitnan')';% calculate mean in cold sesons
dVq_cold_std = std(dVq_cold,0,2,'omitnan')';% calcualte std
dVq_cold_median = median(dVq_cold,2,'omitnan')';% calculate median
tf_nan= isnan(dVq_cold_mean);
p_dVq_cold = (Vq - Vq2)./((Vq+Vq2)./2).*100;% % difference
% p_dVq_cold = (Vq - Vq2)./Vq2.*100;
p_dVq_cold_mean = mean(p_dVq_cold,2,'omitnan')';
p_dVq_cold_std = std(p_dVq_cold,0,2,'omitnan')';
p_dVq_cold_median = median(p_dVq_cold,2,'omitnan')';


%% general difference
figure;
subplot(1,2,1);hold all;
boundedline(dVq_mean(:,~tf_nan),h_interp(:,~tf_nan), dVq_std(:,~tf_nan),'alpha','orientation', 'horiz');% 
plot(dVq_median, h_interp, 'k--');
legend({'1\sigma','mean','median'});
xlabel('GEM-MACH - MAX-DOAS [ppbv]');
ylabel('Altitude [km]');
title('NO_2 profile difference');
subplot(1,2,2);
boundedline(p_dVq_mean(:,~tf_nan),h_interp(:,~tf_nan), p_dVq_std(:,~tf_nan),'alpha','orientation', 'horiz');% 
plot(p_dVq_median, h_interp, 'k--');
legend({'1\sigma','mean','median'});
xlabel('GEM-MACH - MAX-DOAS [%]');
ylabel('Altitude [km]');
title('NO_2 profile rel. difference');

fig_name = ['Profile_vertical_difference_' target];
print_setting(1/2,save_fig,[plot_path fig_name]);

%% warm/cold seasons difference
figure;
c_map = lines(2);
subplot(2,2,1);hold all;
boundedline(dVq_warm_mean(:,~tf_nan),h_interp(:,~tf_nan), dVq_warm_std(:,~tf_nan),'alpha','orientation', 'horiz','cmap', c_map(2,:));% 
plot(dVq_warm_median, h_interp, 'k--');
legend({'1\sigma','mean','median'});
xlabel('GEM-MACH - MAX-DOAS [ppbv]');
ylabel('Altitude [km]');
title('NO_2 profile difference [warm seasons]');grid on;
xlim([-6 6]);
subplot(2,2,2);
boundedline(p_dVq_warm_mean(:,~tf_nan),h_interp(:,~tf_nan), p_dVq_warm_std(:,~tf_nan),'alpha','orientation', 'horiz','cmap', c_map(2,:));% 
plot(p_dVq_warm_median, h_interp, 'k--');
legend({'1\sigma','mean','median'});
xlabel('GEM-MACH - MAX-DOAS [%]');
ylabel('Altitude [km]');
title('NO_2 profile rel. difference [warm seasons]');grid on;
xlim([-250 250]);
subplot(2,2,3);hold all;
boundedline(dVq_cold_mean(:,~tf_nan),h_interp(:,~tf_nan), dVq_cold_std(:,~tf_nan),'alpha','orientation', 'horiz','cmap', c_map(1,:));% 
plot(dVq_cold_median, h_interp, 'k--');
legend({'1\sigma','mean','median'});
xlabel('GEM-MACH - MAX-DOAS [ppbv]');
ylabel('Altitude [km]');
title('NO_2 profile difference [cold seasons]');grid on;
xlim([-6 6]);
subplot(2,2,4);
boundedline(p_dVq_cold_mean(:,~tf_nan),h_interp(:,~tf_nan), p_dVq_cold_std(:,~tf_nan),'alpha','orientation', 'horiz','cmap', c_map(1,:));% 
plot(p_dVq_cold_median, h_interp, 'k--');
legend({'1\sigma','mean','median'});
xlabel('GEM-MACH - MAX-DOAS [%]');
ylabel('Altitude [km]');
title('NO_2 profile rel. difference [cold seasons]');grid on;
xlim([-250 250]);

fig_name = ['Profile_vertical_difference_seasons_' target];
print_setting(1/2,save_fig,[plot_path fig_name]);

%% vertical profile difference
figure;
mycolormap = customcolormap(linspace(0,1,11), {'#68011d','#b5172f','#d75f4e','#f7a580','#fedbc9','#f5f9f3','#d5e2f0','#93c5dc','#4295c1','#2265ad','#062e61'});
subplot(2,1,1);
plot_profile(data.UTC,h_interp,dVq,mycolormap);
% ylim([0 1.5]);
caxis([-20 20]);
colorbar;title('GEM-MACH - MAX-DOAS (NO_2 profile difference) [ppbv]');
ylabel('Altitude [km]');
datetick('x','yyyy-mm','keeplimits');
subplot(2,1,2);
plot_profile(data.UTC,h_interp,p_dVq,mycolormap);
% ylim([0 1.5]);
caxis([-200 200]);
colorbar;title('GEM-MACH - MAX-DOAS (NO_2 profile rel. difference) [%]');
ylabel('Altitude [km]');
datetick('x','yyyy-mm','keeplimits');

fig_name = ['Profile_timeserise_difference_' target];
print_setting(1/2,save_fig,[plot_path fig_name]);

%% scatter plots for different heights
figure;hold all;
% height_bins = [0.2,0.4,0.8,1,1.5,2,3,4];
height_bins = [0.5,1,1.5,2,3,4];
c_map = lines(numel(height_bins));
for j = 1:numel(height_bins)
    if j == 1
        tf_hights = (h_interp < height_bins(j));
    else
        tf_hights = (h_interp >= height_bins(j-1)) & (h_interp < height_bins(j));
    end
    Vq_sub = Vq(tf_hights,:);
    Vq2_sub = Vq2(tf_hights,:);
    Vq_sub = reshape(Vq_sub,[],1);
    Vq2_sub = reshape(Vq2_sub,[],1);

    x = Vq_sub;
    y = Vq2_sub;
    tf_xy_nan = isnan(x) | isnan(y);
    x = x(~tf_xy_nan,:);
    y = y(~tf_xy_nan,:);
    if ~(isempty(x) | isempty(y))
%         line_fits_local(x,y);
        fit_color = c_map(j,:);
        [intercept,slope,slope_nlm,mdl_lm,mdl_nlm,N,R] = line_fits_local(x,y,fit_color);
        intercepts(j) = intercept;
        slopes(j) = slope;
        slope_nlms(j) = slope_nlm;
        Ns(j) = N;
        Rs(j) = R;
    end
end
xlabel('GEM_MACH NO_2 profile [vmr]');
ylabel('MAX-DOAS NO_2 profile [vmr]');
fig_name = ['Scatter_for_layers_' target];
print_setting(1/4,save_fig,[plot_path fig_name]);


figure;
xticklabels_1 = {'0-0.5','0.5-1','1-1.5','1.5-2','2-3','3-4'};
subplot(2,2,1);hold all;
plot(slopes);
plot(slope_nlms);
legend({'Simple linear fits','Zero-intercept fits'},'Location','southeast');
xticks(1:numel(height_bins)); xticklabels(xticklabels_1);
ylabel('Slope');
ylim([-1 2]);
xlabel('Altitude range [km]');grid on;

subplot(2,2,2);
plot(intercepts);
legend({'Simple linear fits'});
xticks(1:numel(height_bins)); xticklabels(xticklabels_1);
ylabel('Intercept');
ylim([0 2]);
xlabel('Altitude range [km]');grid on;

subplot(2,2,3);
plot(Ns,'k--');
xticks(1:numel(height_bins)); xticklabels(xticklabels_1);
ylabel('N');
% ylim([0 100]);
xlabel('Altitude range [km]');grid on;


subplot(2,2,4);
plot(Rs,'k--');
xticks(1:numel(height_bins)); xticklabels(xticklabels_1);
ylabel('R');
ylim([-0.5 1]);
xlabel('Altitude range [km]');grid on;

fig_name = ['Scatter_for_layers_summary_' target];
print_setting(1/2,save_fig,[plot_path fig_name]);

%% diurnal difference

N_grid = size(Vq);
N_grid = N_grid(1);% no. of height grids
data.LST = data.UTC -hours(5);

figure;
cmap =  flipud(jet(N_grid-1));
subplot(2,2,1);
for i = 2:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq(i,tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-2 14]);grid on;
title('GEM-MACH warm seasons');
ylabel('NO_2 ppbv');
xlabel('LST Hour');


subplot(2,2,3);hold all;
for i = 2:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq2(i,tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-2 14]);grid on;
title('MAX-DOAS warm seasons');
ylabel('NO_2 ppbv');
xlabel('LST Hour');

subplot(2,2,2);
for i = 2:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(~tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq(i,~tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-2 14]);grid on;
title('GEM-MACH cold seasons');
ylabel('NO_2 ppbv');
xlabel('LST Hour');


subplot(2,2,4);hold all;
for i = 2:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(~tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq2(i,~tf_warm)';
    table_1layer.LST = [];
    tf_0hr = table_1layer.Hour == 0;
    table_1layer(tf_0hr,:) = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-2 14]);grid on;
title('MAX-DOAS cold seasons');
ylabel('NO_2 ppbv');
xlabel('LST Hour');

fig_name = ['Diurnal_profile_layers_' target];
print_setting(1/2,save_fig,[plot_path fig_name]);

%% >= 1 km
N_grid = size(Vq);
N_grid = N_grid(1);% no. of height grids
data.LST = data.UTC -hours(5);

figure;
cmap =  flipud(jet(N_grid-1));
subplot(2,2,1);
for i = 6:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq(i,tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-2 4]);grid on;
title('GEM-MACH warm seasons [>= 1km]');
ylabel('NO_2 ppbv');
xlabel('LST Hour');


subplot(2,2,3);hold all;
for i = 6:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq2(i,tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-2 4]);grid on;
title('MAX-DOAS warm seasons [>= 1km]');
ylabel('NO_2 ppbv');
xlabel('LST Hour');

subplot(2,2,2);
for i = 6:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(~tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq(i,~tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-2 4]);grid on;
title('GEM-MACH cold seasons [>= 1km]');
ylabel('NO_2 ppbv');
xlabel('LST Hour');


subplot(2,2,4);hold all;
for i = 6:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(~tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq2(i,~tf_warm)';
    table_1layer.LST = [];
    tf_0hr = table_1layer.Hour == 0;
    table_1layer(tf_0hr,:) = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-2 4]);grid on;
title('MAX-DOAS cold seasons [>= 1km]');
ylabel('NO_2 ppbv');
xlabel('LST Hour');

fig_name = ['Diurnal_profile_layers_1km_and_up_' target];
print_setting(1/2,save_fig,[plot_path fig_name]);

%% >= 2 km
N_grid = size(Vq);
N_grid = N_grid(1);% no. of height grids
data.LST = data.UTC -hours(5);

figure;
cmap =  flipud(jet(N_grid-1));
subplot(2,2,1);
for i = 11:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq(i,tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-0.5 1.2]);grid on;
title('GEM-MACH warm seasons [>= 2km]');
ylabel('NO_2 ppbv');
xlabel('LST Hour');


subplot(2,2,3);hold all;
for i = 11:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq2(i,tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-0.5 1.2]);grid on;
title('MAX-DOAS warm seasons [>= 2km]');
ylabel('NO_2 ppbv');
xlabel('LST Hour');

subplot(2,2,2);
for i = 11:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(~tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq(i,~tf_warm)';
    table_1layer.LST = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-0.5 1.2]);grid on;
title('GEM-MACH cold seasons [>= 2km]');
ylabel('NO_2 ppbv');
xlabel('LST Hour');


subplot(2,2,4);hold all;
for i = 11:N_grid-1
    table_1layer = table;
    table_1layer.LST = data.LST(~tf_warm,:);
    table_1layer.Hour = table_1layer.LST.Hour;
    table_1layer.Vq_1layer = Vq2(i,~tf_warm)';
    table_1layer.LST = [];
    tf_0hr = table_1layer.Hour == 0;
    table_1layer(tf_0hr,:) = [];
    hourly_mean = grpstats(table_1layer,'Hour','mean');
    hourly_std = grpstats(table_1layer,'Hour','std');
    boundedline(hourly_mean.Hour, hourly_mean.mean_Vq_1layer, hourly_std.std_Vq_1layer ,'alpha','orientation', 'vert','cmap', cmap(i,:));% 
end
xlim([0 23]);ylim([-0.5 1.2]);grid on;
title('MAX-DOAS cold seasons [>= 2km]');
ylabel('NO_2 ppbv');
xlabel('LST Hour');

fig_name = ['Diurnal_profile_layers_2km_and_up_' target];
print_setting(1/2,save_fig,[plot_path fig_name]);

function plot_profile(UTC,h,profile,mycolormap)
% mycolormap = customcolormap(linspace(0,1,12), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691','#ffffff',});
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

function [intercept,slope,slope_nlm,mdl_lm,mdl_nlm,N,R] = line_fits_local(x,y,fit_color)
% simple linear fit

if nargin ==2
    only_simple = false;
end

mdl_lm = fitlm(x,y);

myfun = @(b,x)b(1)*x;
beta =[1];
mdl_nlm = fitnlm(x,y,myfun,beta);
% mdl_nlm = fitlm(x,y,'Intercept',false);% this is another way to fit with intercept set to 0
% figure;hold all;
plot(x,y,'.','Color',fit_color);
plot(x,predict(mdl_lm,x),'color',fit_color);
% plot(x,predict(mdl_nlm,x),'color',fit_color);

% plot(x,x,'k--');
intercept = mdl_lm.Coefficients.Estimate(1);
slope = mdl_lm.Coefficients.Estimate(2);
slope_nlm = mdl_nlm.Coefficients.Estimate(1);
textbp(['y = ' num2str(slope) '*x + ' num2str(intercept)],'color',fit_color);

%     textbp(['y = ' num2str(slope_nlm) '*x'],'color',fit_color);

% textbp(['N = ' num2str(numel(x)) ]);
% textbp(['R = ' num2str(corr(x,y)) ]);
N = numel(x);
R = corr(x,y);

function [intercept,slope,slope_nlm,mdl_lm,mdl_nlm] = line_fits_local2(x,y,only_simple)
% simple linear fit

if nargin ==2
    only_simple = false;
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
plot(x,x,'k--');
intercept = mdl_lm.Coefficients.Estimate(1);
slope = mdl_lm.Coefficients.Estimate(2);
slope_nlm = mdl_nlm.Coefficients.Estimate(1);
textbp(['y = ' num2str(slope) '*x + ' num2str(intercept)],'color',[0.1 0.1 0.9]);
if ~only_simple
    textbp(['y = ' num2str(slope_nlm) '*x'],'color',[0.9 0.1 0.5]);
end
textbp(['N = ' num2str(numel(x)) ]);
textbp(['R = ' num2str(corr(x,y)) ]);
