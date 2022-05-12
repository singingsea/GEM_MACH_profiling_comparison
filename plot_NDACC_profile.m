function plot_NDACC_profile()
save_fig = 1;
output_path = 'D:\Projects\EnviMes\output_profiles\';
% instrument = '1694_2'; tragets = {'NO2_428NM','O4_360NM','HCHO_343NM'};tragets_gas = {'NO_2','Aerosol Extinction','HCHO'};
instrument = '1694_1'; tragets = {'NO2_460NM','O4_477NM'};tragets_gas = {'NO_2','Aerosol Extinction'};


years = 2016:2021;

plot_path = [output_path 'plot\' instrument '\'];
mkdir(output_path);mkdir(plot_path);
instrument_label = strrep(instrument,'_','-');

for i = 1:numel(years)% loop over years
    year = num2str(years(i));
    for j = 1:numel(tragets)% loop over different gas
        target = char(tragets(j));
        target_label = strrep(target,'_',' ');% this is for title
        traget_gas = char((tragets_gas(j)));% this is for color bar title
        disp(['Printting: EnviMes' instrument ' ' target ' ' year ]);
        data = [];c_max=[];
        load(['D:\Projects\EnviMes\output_profiles\' instrument '\' target '_' instrument '_' year '.mat']);
        
        % qa
        tf = data.qa >= 2;
        data(tf,:) = [];
        
        if ~isempty(data)
            % interp
            interp_timestamp = min(data.UTC):minutes(15):max(data.UTC);
            interp_method = 'mean';
%             interp_method = 'linear';
            data = table2timetable(data);
            data = retime(data,interp_timestamp,interp_method);
            data.UTC_datenum = datenum(data.UTC);
            h_interp = 0:1:4000;
            X = data.UTC_datenum';
            Y = h';
            V = data.profile';% number density profile [molec/cm3]
            V2 = data.vmr';% vmr profile [ppbv]
            Xq= X;
            Yq = h_interp;
            Vq = interp1(Y,V,Yq);
            Vq2 = interp1(Y,V2,Yq);
            data.LST = data.UTC - hours(5);

            % plotting
            switch traget_gas
                case 'NO_2'
                    c_max = prctile(prctile(Vq,97),97);
%                     c_max2 = prctile(prctile(Vq2,97),97);                    
                    c_max2 = 20;
                case 'HCHO'
                    c_max = prctile(prctile(Vq,99),99);
                    c_max2 = prctile(prctile(Vq2,99),99);
                case 'Aerosol Extinction'
                    c_max = prctile(prctile(Vq,97),97);
                    c_max2 = prctile(prctile(Vq2,97),97);
            end
            
            % profile figure [number density plot]
            figure;hold all;m = 1;            
            %for k = 1:3:12% loop over months, we plot the data quarterly
            for k = 1:1:12% loop over months, we plot the data monthly
                subplot(12,1,m);hold all;m = m +1;
                if k == 1% only print title for first subplot
                    title(['EnviMes' instrument_label ' ' target_label ' profile [' year ']']);
%                      print_setting(1,0,[plot_path ['1']]);
                end
%                 tf = (data.LST.Month >= k) & (data.LST.Month <= k+2);
                tf = (data.LST.Month == k);
                plot_profile(data.LST(tf,:),h_interp,Vq(:,tf), ' ');            
%                 xlim([datenum(datetime(str2num(year),k,1)) datenum(datetime(str2num(year),k+3,1))]);
                xlim([datenum(datetime(str2num(year),k,1)) datenum(datetime(str2num(year),k+1,1))]);
    %         datetick('x','ddd','keeplimits')
                ylim([300 1500]);                
                caxis([0 c_max]);
                
                if k == 6 
                    ylabel('Altitude [m]');
                end
                
                if k == 6
                    hcb = colorbar('Location','manual');
                    hcb.Position = [0.93, 0.1, 0.02, 0.8];%  [left, bottom, width, height]
                    colorTitleHandle = get(hcb,'Title');
                    if contains(target,'O4')
                        titleString = [traget_gas ' [km^-^1]'];
                    else
                        titleString = [traget_gas ' [molec/cm^3]'];
                    end
                    current_position = colorTitleHandle.Position;
                    current_position(1,1) = -10;
                    current_position(1,2) = 550;
                    set(colorTitleHandle ,'String',titleString,'Position',current_position);                            
                end
                
            end            
            fig_name = ['EnviMes' instrument_label '_' target_label '_' year];
            print_setting(1,save_fig,[plot_path fig_name]);
            
            
            % profile figure [vmr plot]
            figure;hold all;m = 1;            
            %for k = 1:3:12% loop over months, we plot the data quarterly
            for k = 1:1:12% loop over months, we plot the data monthly
                subplot(12,1,m);hold all;m = m +1;
                if k == 1% only print title for first subplot
                    title(['EnviMes' instrument_label ' ' target_label ' profile [' year ']']);
%                      print_setting(1,0,[plot_path ['1']]);
                end
%                 tf = (data.LST.Month >= k) & (data.LST.Month <= k+2);
                tf = (data.LST.Month == k);
                plot_profile(data.LST(tf,:),h_interp,Vq2(:,tf), ' ');            
%                 xlim([datenum(datetime(str2num(year),k,1)) datenum(datetime(str2num(year),k+3,1))]);
                xlim([datenum(datetime(str2num(year),k,1)) datenum(datetime(str2num(year),k+1,1))]);
    %         datetick('x','ddd','keeplimits')
                ylim([300 1500]);                
                caxis([0 c_max2]);
                
                if k == 6 
                    ylabel('Altitude [m]');
                end
                
                if k == 6
                    hcb = colorbar('Location','manual');
                    hcb.Position = [0.93, 0.1, 0.02, 0.8];%  [left, bottom, width, height]
                    colorTitleHandle = get(hcb,'Title');
                    if contains(target,'O4')
                        titleString = [traget_gas ' [km^-^1]'];
                    else
                        titleString = [traget_gas ' [ppbv]'];
                    end
                    current_position = colorTitleHandle.Position;
                    current_position(1,1) = -10;
                    current_position(1,2) = 550;
                    set(colorTitleHandle ,'String',titleString,'Position',current_position);                            
                end
                
            end            
            fig_name = ['EnviMes' instrument_label '_' target_label '_' year '_vmr'];
            print_setting(1,save_fig,[plot_path fig_name]);
            
            
            % vcd figure1
            figure;hold all;
            plot(data.UTC,data.tvcd,'.');
            title(['EnviMes' instrument_label ' ' target_label ' trop. column [' year ']']);
            if contains(target,'O4')                
                ylabel(['AOD']);
            else
                ylabel([traget_gas ' trop column [molec/cm^3]']);
            end
            grid on;
            fig_name = ['EnviMes' instrument_label '_' target_label '_' year '_trop_column'];
            print_setting(1/2,save_fig,[plot_path fig_name]);
            
            % vcd figure2 with error bar
            figure;hold all;
            errorbar(datenum(data.UTC),data.tvcd,data.tvcd_error);
            datetick('x','mmm-dd','keeplimits');
            title(['EnviMes' instrument_label ' ' target_label ' trop. column [' year ']']);
            if contains(target,'O4')                
                ylabel(['AOD']);
            else
                ylabel([traget_gas ' trop column [molec/cm^3]']);
            end
            grid on;
            fig_name = ['EnviMes' instrument_label '_' target_label '_' year '_trop_column_with_errorbar'];
            print_setting(1/2,save_fig,[plot_path fig_name]);            
            
        else 
           disp('Warning: empty results! No plots can be made.' ) 
        end
    end
    close all;
end



%%
% figure;hold all;
% plot_profile(data.LST,h,data.profile',['NO_2 [molec/cm^3]']);
% ylim([0 4000]);
% % fig_name = [target '_' instrument '_' year '_raw'];
% % print_setting(1/2,save_fig,[plot_path fig_name]);


% no2 example
% figure;hold all;
% plot_profile(data.LST,h_interp,Vq,['NO_2 [molec/cm^3]']);
% ylim([300 4000]);
% x_min = min(datenum(data.LST));x_max = max(datenum(data.LST));
% xlim([x_min x_max]);
% 
% xlim([datenum('2019-12-14') datenum('2019-12-20 23:00')]);
% datetick('x','ddd','keeplimits')
% ylim([300 1500]);
% caxis([0 0.8e12]);
% title('MAX-DOAS NO_2 profile [2019 Dec. 14 to 20]');
% fig_name = ['MAX-DOAS_NO2_2019_example1'];
% print_setting(1/2,save_fig,[plot_path fig_name]);

% hcho example
% figure;hold all;
% plot_profile(data.LST,h_interp,Vq,['HCHO [molec/cm^3]']);
% ylim([300 4000]);
% x_min = min(datenum(data.LST));x_max = max(datenum(data.LST));
% xlim([x_min x_max]);
% xlim([datenum('2019-12-14') datenum('2019-12-20 23:00')]);
% datetick('x','ddd','keeplimits')
% ylim([300 1500]);
% caxis([0 0.5e11]);
% title('MAX-DOAS HCHO profile [2019 Dec. 14 to 20]');
% fig_name = ['MAX-DOAS_HCHO_2019_example1'];
% print_setting(1/2,save_fig,[plot_path fig_name]);


% % hcho example summer
% figure;hold all;
% plot_profile(data.LST,h_interp,Vq,['HCHO [molec/cm^3]']);
% ylim([300 4000]);
% x_min = min(datenum(data.LST));x_max = max(datenum(data.LST));
% xlim([x_min x_max]);
% xlim([datenum('2019-08-11') datenum('2019-08-17 23:00')]);
% datetick('x','ddd','keeplimits')
% ylim([300 1500]);
% caxis([0 3e11]);
% title('MAX-DOAS HCHO profile [2019 Aug. 11 to 17]');
% fig_name = ['MAX-DOAS_HCHO_2019_example2'];
% print_setting(1/2,save_fig,[plot_path fig_name]);

% no2 example summer
% figure;hold all;
% plot_profile(data.LST,h_interp,Vq,['NO2 [molec/cm^3]']);
% ylim([300 4000]);
% x_min = min(datenum(data.LST));x_max = max(datenum(data.LST));
% xlim([x_min x_max]);
% xlim([datenum('2019-08-11') datenum('2019-08-17 23:00')]);
% datetick('x','ddd','keeplimits')
% ylim([300 1500]);
% caxis([0 0.8e12]);
% title('MAX-DOAS NO_2 profile [2019 Aug. 11 to 17]');
% fig_name = ['MAX-DOAS_NO2_2019_example2'];
% print_setting(1/2,save_fig,[plot_path fig_name]);


%% plot
function plot_profile(UTC,h,profile,colorTitle)
mycolormap = customcolormap(linspace(0,1,12), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691','#ffffff',});
% mycolormap = customcolormap(linspace(0,1,11), {'#a60026','#d83023','#f66e44','#faac5d','#ffdf93','#ffffbd','#def4f9','#abd9e9','#73add2','#4873b5','#313691'});
colormap(mycolormap);


imagesc(datenum(UTC),h,profile);
% contourf(datenum(UTC),h,profile) 
set(gca, 'YDir','normal');
datetick('x','mmm-dd','keeplimits')
% hcb = colorbar;
% colorTitleHandle = get(hcb,'Title');
% titleString = colorTitle;
% set(colorTitleHandle ,'String',titleString);
% ylabel('Altitude [m]');
grid on;