%% from USimagingPhantoms.m

file_dir = '..\data\imagingUS\';


% 180626 optiFibreKnot angled 4
%     file_name = '180626\optifibreKnot_angled4_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_10s39m17h_26-06-18_avg1_2D_raw.SGL';
%     trigger_delay = 5e-6;
%     samples_cut_off = 0;
%     samples_t0_correct = -6;
%     c0 = 1484;

% 180626 polymerLeaf curved
%     file_name = '180626\polymerLeaf2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_51s28m20h_26-06-18_avg1_2D_raw.SGL';
%     trigger_delay = 4e-6;
%     samples_cut_off = 0;
%     samples_t0_correct = -6;
%     c0 = 1484;

% 180629 gel wax tmm
%     file_name = '180629\gelwaxLayers_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_39s51m16h_29-06-18_avg1_2D_raw.SGL';
%     trigger_delay = 0;
%     samples_cut_off = 10;
%     samples_t0_correct = -6;
%     c0 = 1470;
%         % clinical scanner:
%         file_name_cli = '180629\clinical scanner\29-06-2018_17-41-25 DICOM +2 higher gain\17-38-33.dcm';

% 181204 atmm orgasol 1
%     file_name = '181204/atmm_orgasol1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
%     trigger_delay = 0;
%     samples_cut_off = 50;
%     samples_t0_correct = -6;
%     c0 = 1544;

% 180828 pork belly 3
    file_name = '180828\porkBelly3_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_26s20m19h_28-08-18_avg1_2D_raw.SGL';
    trigger_delay = 0;
    samples_cut_off = 10;
    samples_t0_correct = -4;
    c0 = 1460;

% 190114 lymph node (L3)
%     file_name = '190114/lymphNode2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_13s51m16h_14-01-19_avg1_2D_raw.SGL';
%     trigger_delay = 0;
%     samples_cut_off = 50;
%     samples_t0_correct = -6;
%     c0 = 1520;


%% load sensor data

disp(['Loading: ' file_name])

[sensor_data, params] = loadSGL([file_dir file_name]);
% sensor_data = - sensor_data;            % flip for trolley scanner

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;


%% view sensor data

disp(['Viewing: ' file_name])

fig_data = figure;
set(gcf,'Position',[100,50,600,800])
imagesc(squeeze(sensor_data(40,:,:))')
    colormap(gray)
    colorbar
    title(strtok(file_name,'@'),'Interpreter','None')
    xlabel('x axis [dx]')
    ylabel('time [dt]')
    drawnow

% pause

% sensor_data = sensor_data(:,:,1:1000);


%% look at frequency content of data

% [frequency, f_series_avg] = freqSpecSGLavg(sensor_data,1/params.dt);
% 
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


%% run reconstruction

disp(['Reconstructing: ' file_name])

[reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
                            'ZeroPad', 10, ...
                            'Upsample', true, ...
                            'Apodise', false, ...
                            'FreqBandFilter', {}, ... % 10e6, 15e6
                            'FreqLowFilter', {}, ... % 30e6
                            'TimeGainCompensate', {}, ...
                            'EnvelopeDetect', true, ...
                            'LogCompress', 0, ...
                            'SaveImageToFile', false ...
                        );

%% post processing

sliceViewer


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


