
clear all
close all

file_dir = '../data/sensorCharac/';

%% all-hard dielectric FP sensor

% file_name = 'allharddielectric/USTplanar3M5_allharddielectric_avg200@0nm_t0[-1250]_dx[100µm]_dy[100µm]_dt[4ns]_05s19m16h_05-03-18_avg200_2D_raw.SGL';
% dx = 1e-4;
% dy = dx;
% dt = 4e-9;
% t_0 = 5e-6;
% ROI : (40:120,30:110)

%% BF4 FP sensor

% file_name = 'BF4/USTplanar3M5_BF4_avg40@0nm_t0[-850]_dx[100µm]_dy[100µm]_dt[4ns]_53s24m20h_13-03-18_avg40_2D_raw.SGL';
% dx = 1e-4;
% dy = dx;
% dt = 4e-9;
% t_0 = 3.4e-6;
% ROI : (30:120,30:120)

%% new all-hard dielectric FP sensor

file_name = 'allharddielectric-R1/USTplanar3M5_AHD1_avg100@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_40s19m12h_03-10-18_avg100_2D_raw.SGL';
t_0 = 0;


%%
slice_x = 20;
slice_y = 20;

[dataSGL, params] = loadSGL([file_dir file_name]);
% dataSGL = dataSGL(30:120,30:120,:);
viewSGL(dataSGL,params,t_0,slice_x,slice_y)

%%
caxis([0,0.03])

%% get PMAX and TOA in 2d grid and vector

[ PMAX_xy, TOA_xy ] = max(-dataSGL,[],3);
[Nx,Ny] = size(TOA_xy);
x = 1:Nx;
y = 1:Ny;
[X,Y] = meshgrid(x,y);X=X';Y=Y';
Xv = reshape(X,[Nx*Ny 1]);
Yv = reshape(Y,[Nx*Ny 1]);
PMAX_vec = reshape(PMAX_xy,[Nx*Ny 1]);
TOA_vec  = reshape(TOA_xy, [Nx*Ny 1]);

%% uniformity (histogram plot of PMAX)

% figure
hold on
histogram(PMAX_vec,'BinWidth',0.02)
    title('histogram of peak pressure across scanning area')
    xlabel('peak signal amplitude')
    ylabel(['count (out of ' num2str(Nx*Ny) ')'])
    %legend('whole scanning area','within ROI')


%% sensor frequency response BK31

%reference measurement with AHD1 / freq content of generation:
file_dir_AHD1 = 'D:\PROJECT\data\sensorCharac\allharddielectric-R1\191206\';
file_names_AHD1 = {
                'laserGenSPdiffSource_AHD1_SPI_x0_y0_avg1000.txt';
                'laserGenSPdiffSource_AHD1_SPI_x1_y0_avg1000.txt';
                'laserGenSPdiffSource_AHD1_SPI_x0_y1_avg1000.txt';
                'laserGenSPdiffSource_AHD1_SPI_x1_y1_avg1000.txt'
              };

%measurement of BK31 / sensor freq response * freq content of generation:
file_dir_BK31 = 'D:\PROJECT\data\sensorCharac\BK31\191209\';
file_names_BK31 = {
                'laserGenSPdiffSource_BK31_SPI_x0_y0_avg100.txt';
                'laserGenSPdiffSource_BK31_SPI_x1_y0_avg100.txt';
                'laserGenSPdiffSource_BK31_SPI_x0_y1_avg100.txt';
                'laserGenSPdiffSource_BK31_SPI_x1_y1_avg100.txt'
              };

t_0 = 3e-6;
dt = 4e-9;

% time = viewSGLsingle(file_dir,file_name,t_0,'removeDC',true);

t_min = 1;
t_max = 500;

f_series_AHD1_ar = zeros(4,251);
f_series_BK31_ar = zeros(4,251);

sensor_freq_response_ar = zeros(16,251);

for idx_AHD1 = 1:length(file_names_AHD1)

    file_name_AHD1 = file_names_AHD1{idx_AHD1};
    [frequency_AHD1, f_series_AHD1] = freqSpecSGLsingle(file_dir_AHD1,file_name_AHD1,1/dt,t_min,t_max,'NormTime',false,'NormFreq',true,'LineColour','g');
    f_series_AHD1_ar(idx_AHD1,:) = f_series_AHD1;
    
end

for idx_BK31 = 1:length(file_names_BK31)

    file_name_BK31 = file_names_BK31{idx_BK31};
    [frequency_BK31, f_series_BK31] = freqSpecSGLsingle(file_dir_BK31,file_name_BK31,1/dt,t_min,t_max,'NormTime',false,'NormFreq',true,'LineColour','r');
    f_series_BK31_ar(idx_BK31,:) = f_series_BK31;

    assert(isequal(frequency_AHD1, frequency_BK31))

end

idx_sfr = 0;
for idx_AHD1 = 1:length(file_names_AHD1)
    for idx_BK31 = 1:length(file_names_BK31)
%         idx_BK31 = idx_AHD1;
        
        sensor_freq_response = f_series_BK31_ar(idx_BK31,:) ./ f_series_AHD1_ar(idx_AHD1,:);
        
        idx_sfr = idx_sfr + 1;
        sensor_freq_response_ar(idx_sfr,:) = sensor_freq_response;
        
    end
end

% normalisation of each array to common norm
% f_series_AHD1_ar = f_series_AHD1_ar / max(max(f_series_AHD1_ar));
% f_series_BK31_ar = f_series_BK31_ar / max(max(f_series_BK31_ar));
sensor_freq_response_ar = sensor_freq_response_ar / max(max(sensor_freq_response_ar(:,1:100)));

% mean and std for each array
f_series_AHD1_mean = mean(f_series_AHD1_ar,1);
f_series_BK31_mean = mean(f_series_BK31_ar,1);
sensor_freq_response_mean = mean(sensor_freq_response_ar,1);

f_series_AHD1_std = std(f_series_AHD1_ar,[],1);
f_series_BK31_std = std(f_series_BK31_ar,[],1);
sensor_freq_response_std = std(sensor_freq_response_ar,[],1);

% plot normalised data
figure(1002)
set(gcf,'Position',[100 20 700 400])
set(gca,'YScale','log')
hold on
for idx_AHD1 = 1:size(f_series_AHD1_ar,1)
    semilogy(frequency_AHD1/1e6, f_series_AHD1_ar(idx_AHD1,:),'g.','MarkerSize',4)
end
for idx_BK31 = 1:size(f_series_BK31_ar,1)
    semilogy(frequency_BK31/1e6, f_series_BK31_ar(idx_BK31,:),'r.','MarkerSize',4)
end
for idx_sfr = 1:size(sensor_freq_response_ar,1)
    semilogy(frequency_BK31/1e6, sensor_freq_response_ar(idx_sfr,:),'b.','MarkerSize',4)
end
ylim([1e-3 1.6])
xlim([0,80])
xlabel('frequency / MHz')
ylabel('signal amplitude / V')

% plot means & std over data
figure(1002)
hold on
plot(frequency_AHD1/1e6,f_series_AHD1_mean,'k-')
plot(frequency_AHD1/1e6,f_series_AHD1_mean-f_series_AHD1_std,'k--')
plot(frequency_AHD1/1e6,f_series_AHD1_mean+f_series_AHD1_std,'k--')
plot(frequency_BK31/1e6,f_series_BK31_mean,'k-')
plot(frequency_BK31/1e6,f_series_BK31_mean-f_series_BK31_std,'k--')
plot(frequency_BK31/1e6,f_series_BK31_mean+f_series_BK31_std,'k--')

plot(frequency_BK31/1e6,sensor_freq_response_mean,'k-')
plot(frequency_BK31/1e6,sensor_freq_response_mean-sensor_freq_response_std,'k--')
plot(frequency_BK31/1e6,sensor_freq_response_mean+sensor_freq_response_std,'k--')

h = zeros(3, 1);
h(1) = plot(NaN,NaN,'g:');
h(2) = plot(NaN,NaN,'r:');
h(3) = plot(NaN,NaN,'b:');
legend(h, 'AHD1 reference','BK31 sensor','BK31 / AHD1');
title('sensor frequency response (4*4 = 16 avg)')
% title('sensor frequency response (4*1to1 = 4 avg)')
set(gca,'FontSize',13)


%% plot normalised sensor freq response

figure
set(gcf,'Position',[100 20 700 400])
plot(frequency_BK31/1e6,sensor_freq_response_mean/max(sensor_freq_response_mean(1:100)),'k-'), ylim([1e-3 1.6])
    set(gca,'YScale','log')
    xlim([0,80])
    ylim([0.04,1.2])
    xlabel('frequency / MHz')
    ylabel('signal amplitude / V')
    hold on


%% save normalised sensor freq response

frequency = frequency_BK31;
sensor_freq_response_mean_norm = sensor_freq_response_mean/max(sensor_freq_response_mean(1:100));

save('D:\PROJECT\data\sensorCharac\BK31\191209\sensor_frequency_response.mat', 'frequency', 'sensor_freq_response_mean_norm')


%% search for equivalent gaussian bandpass filter for simulation

centre_freq = 1e6;      % [Hz]
bandwidth   = 66e6;     % [Hz]
f = frequency_BK31;     % [Hz]

% copied from gaussianFilter (which is used for freq bandpass filtering in simul_scattTMM_analysis
mean = centre_freq;
variance = (bandwidth / (2 * sqrt(2 * log(2)))).^2;
magnitude = 1;
gauss_filter = max(gaussian(f, magnitude, mean, variance), gaussian(f, magnitude, -mean, variance));


gcf
hold on
plot(f/1e6,gauss_filter)




