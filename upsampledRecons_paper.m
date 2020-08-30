%% from USimagingPhantoms.m

file_dir = '..\data\imagingUS\';

% 191126 resolution27umPlanar (to look at sensor data only . reconstruction in separate script)
%     file_name = '191126\resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_x-05_z0@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_29s48m13h_26-11-19_avg1_1D_raw.SGL'; trigger_delay = 0;
%     file_name = '191126\resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_x-05_z3@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_13s32m14h_26-11-19_avg1_1D_raw.SGL'; trigger_delay = 3e-6;
%     file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\USimaging\191126 resolution27umPlanar BK31[CNT] trolley scrambled fibre centralised parallel phantom\';

% 180626 optiFibreKnot angled 4
    file_name = '180626\optifibreKnot_angled4_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_10s39m17h_26-06-18_avg1_2D_raw.SGL';
    trigger_delay = 5e-6;
    samples_cut_off = 0;
    samples_t0_correct = -6;
    c0 = 1484;
    trim_tz = 1:700;
    file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\USimaging\180626 optifibreKnot angled BK31[CNT]\';

% 180626 polymerLeaf curved
    file_name = '180626\polymerLeaf2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_51s28m20h_26-06-18_avg1_2D_raw.SGL';
    trigger_delay = 4e-6;
    samples_cut_off = 0;
    samples_t0_correct = -6;
    c0 = 1484;
    trim_tz = 1:800;
    file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\USimaging\180626 polymerLeaf BK31[CNT]\';

% 180629 gel wax tmm
    file_name = '180629\gelwaxLayers_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_39s51m16h_29-06-18_avg1_2D_raw.SGL';
    trigger_delay = 0;
    samples_cut_off = 10;
    samples_t0_correct = -6;
    c0 = 1470;
    trim_tz = 100:1200; % 100:2200;
    file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\USimaging\180629 gelwaxLayers BK31[CNT]\';
%         % clinical scanner:
%         file_name_cli = '180629\clinical scanner\29-06-2018_17-41-25 DICOM +2 higher gain\17-38-33.dcm';

% 181204 atmm orgasol 1
    file_name = '181204/atmm_orgasol1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
    trigger_delay = 0;
    samples_cut_off = 50;
    samples_t0_correct = -6;
    c0 = 1544;
    trim_tz = 200:2300;
    file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\USimaging\181204 atmm orgasol BK31[CNT]\';

% 180828 pork belly 3
    file_name = '180828\porkBelly3_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_26s20m19h_28-08-18_avg1_2D_raw.SGL';
    trigger_delay = 0;
    samples_cut_off = 10;
    samples_t0_correct = -4;
    c0 = 1460;
    trim_tz = 150:550;
    file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\USimaging\180828 porkBelly3 BK31[CNT]\';

% 190114 lymph node (L3)
    file_name = '190114/lymphNode2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_13s51m16h_14-01-19_avg1_2D_raw.SGL';
    trigger_delay = 0;
    samples_cut_off = 50;
    samples_t0_correct = -6;
    c0 = 1520;
    trim_tz = 700:1500;
    file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\USimaging\190114 lymph node 2 BK31[CNT]\';


%% load sensor data

disp(['Loading: ' file_name])

[sensor_data, params] = loadSGL([file_dir file_name]);
% sensor_data = - sensor_data;            % flip for trolley scanner

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;


%% view sensor data

kgrid = kWaveGrid(params.Nx, params.dx, params.Ny, params.dy);
t_array = trigger_delay + linspace(1,params.Nt,params.Nt)*params.dt;

% trim_tz = 1:1000;

sensor_data = sensor_data(:,:,trim_tz);
t_array = t_array(trim_tz);

disp(['Viewing: ' file_name])

half_x = round(params.Nx/2);
half_y = round(params.Ny/2);

fig_data = figure;
set(gcf,'Position',[100,50,600,800])
imagesc(kgrid.x_vec*1e3, t_array*1e6, squeeze(sensor_data(half_x,:,:))')
    colormap(gray)
    colorbar
    title(strtok(file_name,'@'),'Interpreter','None')
    xlabel('x axis [mm]')
    ylabel('time [\mum]')
    set(gca,'FontSize',13)
    drawnow

fig_data_1d = figure;
set(gcf,'Position',[450,100,800,500])
plot(t_array*1e6, squeeze(sensor_data(half_x,half_y,:)))
    title(strtok(file_name,'@'),'Interpreter','None')
    xlabel('time [\mum]')
    ylabel('signal amplitude [V]')
    set(gca,'FontSize',13)
    ylim([-0.1,0.14])
    drawnow

% pause


%% look at frequency content of data

[frequency, f_series_avg] = freqSpecSGLavg(sensor_data,1/params.dt);

% legend('optiFibreKnot','polymerLeaf','gelwaxLayers','agarTMMorgasol','porkBelly','lymphNode')

% file_name_fig = [file_dir_figs strtok(file_name(8:end),'@')];
% saveas(gcf, [file_name_fig '_data_freq_avg.fig'])
% saveas(gcf, [file_name_fig '_data_freq_avg.jpg'])

% ---

% f_series_avg = zeros(1,params.Nt);
%  
% figure
%     hold on
%     title('single frequency spectrum')
%     xlabel('frequency / MHz')
%     ylabel('signal amplitude / V')
%         
% for slice_x = round(params.Nx/2):params.Nx
%     for slice_y = round(params.Ny/2):params.Ny
%         t_series = squeeze(sensor_data(slice_x,slice_y,:));
%         [frequency, f_series] = spect(t_series,1/params.dt); %,'Window','Cosine');
%         f_series_avg = f_series_avg + f_series;
%         
%         semilogy(frequency/1e6, f_series/max(f_series))
%         drawnow
% 
%     end
% end
% 
% f_series_avg = f_series_avg / (Nx*Ny);


%% 2dfft sensor data
% 
% sensor.data   = squeeze(sensor_data(half_x,:,:));
% % sensor.data   = squeeze(mean(sensor_data,1));
% sensor.params = params;
% 
% [size_x, size_t] = size(sensor.data);
% 
% size_t_fft = round((size_t+1)/2);
% size_x_fft = round((size_x+1)/2);
% 
% sensor_data_fftT  = zeros(size_x    , size_t_fft);
% sensor_data_fftTX = zeros(size_x_fft, size_t_fft);
% 
% for x = 1 : sensor.params.Ny
%     [freqT, sensor_data_fftT(x,:) ] = spect(sensor.data(x,:), 1/sensor.params.dt);
% end
% 
% for t = 1 : size_t_fft
%     [freqX, sensor_data_fftTX(:,t)] = spect(sensor_data_fftT(:,t), 1/sensor.params.dx);
% end
% 
% % % cut out super high freq
% % omegarange = 1:round(length(freqT)/2);
% % freqT = freqT(omegarange);
% % sensor_data_fftTX = sensor_data_fftTX(:,omegarange);
% 
% x_min = 4;
% fig_2dfft = figure('Position',[300,300,750,450]);
% imagesc(freqT/1e6, freqX(x_min:end)/1e3, sensor_data_fftTX(x_min:end,:))
%     title(strtok(file_name,'@'),'Interpreter','None')
%     xlabel('Temporal frequency \omega [MHz]')
%     ylabel('Spatial frequency k_x [mm^{-1}]')
% %     xlim([0,70])
% %     ylim([0,20])
%     colorbar
%     set(gca,'FontSize',13)
%     % caxis([-0.5e-5, 1e-5])
% %     caxis([0,3e-5])
% %     caxis([0,1e-5])


%% save figs
% 
% file_name_fig = [file_dir_figs strtok(file_name(8:end),'@')];
% 
% saveas(fig_data   , [file_name_fig '_data.fig'])
% saveas(fig_data   , [file_name_fig '_data.jpg'])
% saveas(fig_data_1d, [file_name_fig '_data_1d.fig'])
% saveas(fig_data_1d, [file_name_fig '_data_1d.jpg'])
% saveas(fig_2dfft  , [file_name_fig '_data_fft.fig'])
% saveas(fig_2dfft  , [file_name_fig '_data_fft.jpg'])


%% run reconstruction
% 
% disp(['Reconstructing: ' file_name])
% 
% [reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
%                             'ZeroPad', 10, ...
%                             'Upsample', true, ...
%                             'Apodise', false, ...
%                             'FreqBandFilter', {}, ... % 10e6, 15e6
%                             'FreqLowFilter', {}, ... % 30e6
%                             'TimeGainCompensate', {}, ...
%                             'EnvelopeDetect', true, ...
%                             'LogCompress', 0, ...
%                             'SaveImageToFile', false ...
%                         );


%% post processing
% 
% sliceViewer


%% log compress images
% 
% file_dir  = 'D:\PROJECT\code\recon_data\USIPAPER backups\processed - trimmed interpolated tgc log\';
% 
% % file_name = 'optifibreKnot_angled4_BK31[CNT]_fc30_trimmed-15-15-3_interpolated-25-25-10';
% % file_name = 'polymerLeaf2_BK31[CNT]_fc30_trimmed-15-15-3.5_interpolated-25-25-10';
% % file_name = 'gelwaxLayers_BK31[CNT]_fc30_trimmed-15-2-15.5_interpolated-50-50-10_tgc-100';
% % file_name = 'atmm_orgasol1_BK31[CNT]_compound_trimmed-10-2-8_interpolated-25-25-10_tgc150-50';
% % file_name = 'porkBelly3_BK31[CNT]_trimmed-15-2-11_interpolated-25-25-25_tgc150';
% % file_name = 'porkBelly3_BK31[CNT]_trimmed-15-2-11_interpolated-25-25-25_tgc200';
% file_name = 'lymphNode2_BK31[CNT]_f10_bw15_trimmed-13-15-4_interpolated-25-25-10_tgc150';
% 
% % load volume data
% load([file_dir file_name '.mat'],'volume_data','volume_spacing')
% 
% % scale volume data between 0 and 1
% volume_data = (volume_data - min(volume_data(:))) / (max(volume_data(:)) - min(volume_data(:)));
% 
% % log compress volume data with factor 20
% volume_data = 20 * log10(volume_data);
% 
% % save log compressed volume data
% save([file_dir file_name '_logcompressed.mat'],'volume_data','volume_spacing')


%% get dB value of range of amplitudes (after 0->1 scaling)

file_dir  = 'D:\PROJECT\code\recon_data\USIPAPER backups\processed - trimmed interpolated tgc log\';

% file_name = 'optifibreKnot_angled4_BK31[CNT]_fc30_trimmed-15-15-3_interpolated-25-25-10'; % max 1, 2, 3
file_name = 'polymerLeaf2_BK31[CNT]_fc30_trimmed-15-15-3.5_interpolated-25-25-10'; % max 1, 2, 3
% file_name = 'gelwaxLayers_BK31[CNT]_fc30_trimmed-15-2-15.5_interpolated-50-50-10_tgc-100'; % mean 2
% file_name = 'atmm_orgasol1_BK31[CNT]_compound_trimmed-10-2-8_interpolated-25-25-10_tgc150-50'; % mean 2
% file_name = 'porkBelly3_BK31[CNT]_trimmed-15-2-11_interpolated-25-25-25_tgc200'; % mean 2
% file_name = 'lymphNode2_BK31[CNT]_f10_bw15_trimmed-13-15-4_interpolated-25-25-10_tgc150'; % max 1, 2, 3 & slices 120:200 (+40:+40)

% load volume data
load([file_dir file_name '.mat'],'volume_data','volume_spacing')

% scale volume data between 0 and 1
volume_data = (volume_data - min(volume_data(:))) / (max(volume_data(:)) - min(volume_data(:)));

% max or mean projection along dim (see above)
MIP = mean(volume_data,2);
% MIP = max(volume_data,[],2); % (:,320:400,:)

% display min max values
disp([min(volume_data(:)),max(volume_data(:))])
disp([min(MIP(:)),max(MIP(:))])

% dB value
dBvalue = 20*log10(max(MIP(:))/min(MIP(:)));
disp(dBvalue)


%% fly through videos

% % data = load('D:\PROJECT\code\recon_data\porkBelly3_BK31[CNT].mat');
% % data = load('D:\PROJECT\code\recon_data\porkBelly3_BK31[CNT]_trimmed_tgc.mat');
% % data = load('D:\PROJECT\code\recon_data\porkBelly3_BK31[CNT]_trimmed_tgc_interp.mat');
% % data = load('D:\PROJECT\code\recon_data\lymphNode2_BK31[CNT]_f10_bw15_trimmed_tgc.mat');
% data = load('D:\PROJECT\code\recon_data\lymphNode2_BK31[CNT]_f10_bw15_trimmed_tgc_interp.mat');
% 
% image = data.volume_data;
% voxsz = data.volume_spacing;
% 
% % trim in depth & y
% %maxdepth = round( 5e-3/voxsz(3));
% miny     = round( 2e-3/voxsz(2)); % 40;  %
% maxy     = round(12e-3/voxsz(2)); % 240; %
% image    = image(:,miny:maxy,:);    %1:maxdepth);
% 
% % make spatial axes
% x_axis = (0:size(image,1)-1)*voxsz(1)*1e3;    % [mm]
% z_axis = (0:size(image,3)-1)*voxsz(3)*1e3;    % [mm]
% 
% % slice & frame parameters
% slicethickness = round(1e-3/voxsz(2));      % [voxel]
% num_frames     = size(image,2) - slicethickness + 1;
% 
% % open video
% % vidObj = VideoWriter('D:\PROJECT\figures\_Matlab figs\USimaging\180828 porkBelly3 BK31[CNT]\porkBelly3_fly_zx_USIPAPER.avi');   % ,'Uncompressed AVI');
% vidObj = VideoWriter('D:\PROJECT\figures\_Matlab figs\USimaging\190114 lymph node 2 BK31[CNT]\lymphNode2_fly_zx_USIPAPER.avi');   % ,'Uncompressed AVI');
% vidObj.FrameRate = slicethickness;
% open(vidObj);
% 
% % open figure
% figure('Position',[300,300,600,300]) % PORK [300,300,600,450]) % LYMPH
% 
% % add frames to video
% for frame = 1 : num_frames
%     ypos = frame*voxsz(2)*1e3;
%     meanIP = squeeze(mean(image(:,frame:frame+slicethickness-1,:),2));
%     gcf
%     imagesc(x_axis,z_axis,meanIP')
%         axis image
%         colormap(gray)
%         brighten(0.2)
%         xlabel('x axis [mm]')
%         ylabel('depth z [mm]')
%         %annotation('textbox',[0.1,0.1,0.3,0.3],'String',['y = ']
%         text(10.5,3.5,['y = ' sprintf('%0.1f', ypos) ' mm'],'FontSize',14,'FontName','Times New Roman','Color',[1 1 1])     % PORK 11.5,4.5 / 11.5,10.5 % LYMPH 10.5,3.5
%         set(gca,'FontSize',14)
%         set(gca,'FontName','Times New Roman')
%     currFrame = getframe(gcf);
%     writeVideo(vidObj,currFrame);
%     %pause(0.01)
% end
% 
% % check that all frames added
% assert(size(image,2) == frame+slicethickness-1);
% 
% % close video
% close(vidObj);


