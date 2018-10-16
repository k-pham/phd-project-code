% clear all
% 
% file_dir = '../data/heatingEffects';
% 
% %file_name = '170707/USTplanar3M5_laserON.txt';
% %file_name = '170710/USTplanar3M5_laserON_x2_y0.txt';
% %file_name = '170710/USTplanar3M5_laserON_x1_y0.txt';
% %file_name = '170710/USTplanar3M5_laserON_x0_y0.txt';
% %file_name = '170710/USTplanar3M5_laserON_x1_y0_2.txt';
% %file_name = '170710/USTplanar3M5_laserON_x1_y0_3.txt';
% %file_name = '170710/USTplanar3M5_laserON_x1_y0_4.txt';
% 
% % time scale run 60.1 seconds
% %file_name = '170711/timescale_500dt.txt';
% 
% % reach actual SS: leave laser ON for long then OFF for long
% %file_name = '170711/USTplanar3M5_laserON_x0_y0_longSS.txt';
% 
% % measure UST in thermal equilibrium
% %file_name = '170711/USTplanar3M5_laserON_x0_y0_autotrackbias.txt';
% 
% % measure laser gen US reflection in thermal equilibrium
% %file_name = '170717/laserGenUS_waterReflection.txt';
% 
% % long term AC and DC triggered on laser
% %file_name = '170721/laserGenUS_waterReflection_short1.txt';
% %file_name = '170721/laserGenUS_waterReflection_short2.txt';
% %file_name = '170721/laserGenUS_waterReflection_long.txt';
% 
% %file_name = '171003/heatingEffects_BA56_longTerm_trigger[pulser].txt';
% %file_name = '171003/heatingEffects_BA56_pulse2pulse_trigger[laser].txt';
% %file_name = '171003/heatingEffects_BA56_pulse2pulse_trigger[laser]_2.txt';
% 
% %file_name = '171006/heatingEffects_BA56_pulse2pulse_trigger[laser].txt';
% %file_name = '171006/heatingEffects_BA56_pulse2pulse_trigger[laser]_2.txt';
% 
% %file_name = '171027/heatingEffects_BA52[thin]_pulse2pulse_trigger[laser].txt';
% %file_name = '171027/heatingEffects_BA52[thin]_pulse2pulse_trigger[laser]_2.txt';
% 
% %file_name = '171027/heatingEffects_BA56[thick]_pulse2pulse_trigger[laser].txt';
% %file_name = '171027/heatingEffects_BA56[thick]_pulse2pulse_trigger[laser]_2.txt';
% 
% % in BEAM CENTRE !!! and changing PRF of laser:
% %file_name = '171106/heatingEffects_BA56_pulse2pulse_beamcentre_20Hz.txt';
% %file_name = '171106/heatingEffects_BA56_pulse2pulse_beamcentre_20Hz_2.txt';
% %file_name = '171106/heatingEffects_BA56_pulse2pulse_beamcentre_10Hz.txt';
% file_name = '171106/heatingEffects_BA56_pulse2pulse_beamcentre_05Hz.txt';
% 
% 
% file_path = [file_dir '/' file_name];
% 
% %data = importdata(filepath);
% 
% fileID = fopen(file_path,'r');
% formatSpec = '%f %f %f';
% sizeData = [3 Inf];
% data = fscanf(fileID,formatSpec,sizeData);
% fclose(fileID);
% data = data';
% 
% clear file_dir and file_name and file_path
% 
% time = data(:,1);
% vAC  = data(:,2);
% vDC  = data(:,3);

%% hE-long term: extract peak pressure from AC and average DC for every waveform
% 
% num_samples = 500;      % number of samples per acquisition
% dt_acquis = 0.005768;   % dt between acquisitions in seconds
% 
% % number of acquisitions in whole measurement
% num_acquis = int32(length(time) / num_samples);
% 
% % initiate arrays
% ACpeaks = zeros(1,num_acquis);
% DCavgs = zeros(1,num_acquis);
% realtime = zeros(1,num_acquis);
% 
% % for every acquisition find peak AC and average DC
% for i = 1 : num_acquis
%     ACpeak = max(vAC( (i-1)*num_samples+1 : i*num_samples ));
%     DCavg = mean(vDC( (i-1)*num_samples+1 : i*num_samples ));
%     
%     %display( time((i-1)*timesteps+1) )
%     %display( time(i*timesteps) )
%     %display(peak)
%     %display(DCavg)
%     
%     ACpeaks(i) = ACpeak;
%     DCavgs(i) = DCavg;
%     realtime(i) = double(i) * dt_acquis;
% end
% 

%% hE-long term: plot time evolution of peak average DC and AC/pressure vs real time
% 
% figure(1); clf(1)
% plot(realtime,DCavgs,'b')
%     title('vDC average per waveform')
%     xlabel('time [s]')
%     ylabel('average vDC')
%     %xlim([0,190])
%     %ylim([0,2])
% 
% figure(2); clf(2)
% plot(realtime,ACpeaks,'r')
%     title('vAC/pressure peak per waveform')
%     xlabel('time [s]')
%     ylabel('peak vAC')
%     %ylim([0,1])

%% hE-long term: frequency spectrum of DC and AC
% 
% freq_sampling = 1/dt_acquis;
% 
% [frequency, vDC_freq] = spect(DCavgs,freq_sampling);
% [frequency, vAC_freq] = spect(ACpeaks,freq_sampling);       % , 'Plot', [true false]
%     
% 
% figure(5)
% plot(frequency, vDC_freq, 'b')
%     title('frequency component in DC average evolution')
%     xlabel('frequency [Hz]')
%     ylabel('signal amplitude')
%     %xlim([0,40])
%     
% figure(6)
% plot(frequency, vAC_freq, 'r')
%     title('frequency component in AC peak evolution')
%     xlabel('frequency [Hz]')
%     ylabel('signal amplitude')
%     %xlim([0,40])


%% hE-p2p: plot raw DC and AC consecutive acquisitions
% 
% dt_samples = 100e-9;
% num_samples = 3000000;
% realtime = linspace(0,dt_samples*num_samples,num_samples); %in seconds
% realtime = realtime*1e3;    % in ms
% %realtime = realtime - 1;
% vDC = vDC + 0.5;
% 
% figure(6);% clf(6)
% plot(realtime,vDC(1:num_samples),'r')
%     title('vDC in consecutive acquisitions')
%     xlabel('time [ms]')
%     ylabel('vDC')
%     %xlim([-1,9])
% 
% figure(4); clf(4)
% plot(vAC,'r')
%     title('vAC in consecutive acquisitions')
%     xlabel('time arbitrary')
%     ylabel('vAC')


%% compare pre-tuning maps when cold and heated: move scanning area

% xmin = '-5'; xmax = '5'; % middle
% ymin = '-5'; ymax = '5';

% xmin = '-7'; xmax = '3'; % lower x (use this to see beam centre)
% ymin = '-5'; ymax = '5';

% xmin = '-3'; xmax = '7'; % higher x
% ymin = '-5'; ymax = '5';

% xmin = '-5'; xmax = '5'; % lower y
% ymin = '-7'; ymax = '3';

% xmin = '-5'; xmax = '5'; % higher y (use this to see beam centre)
% ymin = '-3'; ymax = '7';

% file_dir = '../data/heatingEffects/';
% file_dir_cold = [ '171106/heatingEffects_BA56_pretuning_cold_x[',xmin,',',xmax,']_y[',ymin,',',ymax,']/' ];
% file_dir_heat = [ '171106/heatingEffects_BA56_pretuning_heat_x[',xmin,',',xmax,']_y[',ymin,',',ymax,']/' ];
% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)

% figure(6)
% subplot(3,3,7)
% caxis([0,1])

%% compare pre-tuning maps when cold and heated: change PRF of exc.laser

% file_dir = '../data/heatingEffects/';
% file_dir_cold = '171106/heatingEffects_BA56_pretuning_x[-7,3]_y[-3,7]_cold/';
% file_dir_heat = '171106/heatingEffects_BA56_pretuning_x[-7,3]_y[-3,7]_heat05Hz/';
% %file_dir_heat = '171106/heatingEffects_BA56_pretuning_x[-7,3]_y[-3,7]_heat10Hz/';
% %file_dir_heat = '171106/heatingEffects_BA56_pretuning_x[-7,3]_y[-3,7]_heat20Hz/';

% repeat with sensor rotated by 45^\circ on 171107
% rotated = '_rotated45/';
% file_dir_cold = [file_dir_cold(1:5) '7' file_dir_cold(7:end-1) rotated];
% file_dir_heat = [file_dir_heat(1:5) '7' file_dir_heat(7:end-1) rotated];

% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)

% figure(6)
% subplot(3,3,7)
% caxis([0,1])

%% laser-pulser-triggering (LPT): compare pulse shap/wave field to pulser-triggering (PT)
% mapped 3.5MHz planar UST with harddielectric sensor for both trigger mechanisms
% 
% file_dir = '../data/heatingEffects/';
% file_name_PT  = '171204/laser_pulser_triggering/USTplanar3M5_harddielectric_trigger[pulser]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_04s25m20h_04-12-17_avg8_2D_raw.SGL';
% file_name_LPT = '171204/laser_pulser_triggering/USTplanar3M5_harddielectric_trigger[laser_pulser]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_14s58m21h_04-12-17_avg8_2D_raw.SGL';
% 
% dataSGL_PT = loadSGL( [file_dir file_name_PT] );
% dataSGL_LPT = loadSGL( [file_dir file_name_LPT] );
% 
% % viewing and comparing single time series
% slice_x = 50;
% slice_y = 50;
% dx = 1e-4;
% dy = dx;
% dt = 20e-9;
% t_0 = 0;
% t_min = 1;
% t_max = size(dataSGL_PT,3);
% 
% viewSGL(dataSGL_PT ,dx,dy,dt,t_0,slice_x,slice_y)
% viewSGL(dataSGL_LPT,dx,dy,dt,t_0,slice_x,slice_y)
% 
% %viewing and comparing single frequency series
% [frequency,f_series_PT]  = freqSpecSGL(dataSGL_PT ,1/dt,slice_x,slice_y,t_min,t_max);
% [frequency,f_series_LPT] = freqSpecSGL(dataSGL_LPT,1/dt,slice_x,slice_y,t_min,t_max);
% plot(frequency/1e6,f_series_LPT-f_series_PT,'k--')
% legend('pulser triggering','laser-pulseGen-pulser triggering','difference')
% 
% % movie comparing time series for USTplanar3M5 mappings using different triggering methods
% Nx = size(dataSGL_PT,1);
% Ny = size(dataSGL_PT,2);
% Nt = size(dataSGL_PT,3);
%     
% kgrid = kWaveGrid(Nx, dx, Ny, dy);
% time = (t_0+linspace(dt,Nt*dt,Nt))*1e6;
% difference = dataSGL_LPT - dataSGL_PT;
% 
% figure
% for slice_x = 1:Nx
%     for slice_y = 1:Ny
%         plot(time,squeeze(dataSGL_PT(slice_x,slice_y,:)),'b')
%         hold on
%         plot(time,squeeze(dataSGL_LPT(slice_x,slice_y,:)),'r')
%         plot(time,squeeze(difference(slice_x,slice_y,:)),'k--')
%             title(['time series at x = ',num2str(kgrid.x_vec(slice_x)*1e3),' mm, y = ',num2str(kgrid.y_vec(slice_y)*1e3),' mm'])
%             xlabel('time / \mu s')
%             ylabel('signal amplitude / V')
%             %legend('pulser triggering','laser-pulseGen-pulser triggering','difference')
%             axis([0,3,-0.2,0.15])
%         hold off
%         drawnow
%         pause(0.001)
%     end
% end

%% pre-tune & map USTplanar3M5 for cold and hot BA52/spraypaint
% 
% clear all 
% close all
% 
% file_dir = '../data/heatingEffects/';
% file_pretune_cold = '171212/pretuning_BA52_cold/';
% file_pretune_heat = '171212/pretuning_BA52_heat/';
% file_acquire_cold  = '171212/USTplanar3M5_BA52_cold@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_38s30m00h_13-12-17_avg4_2D_raw.SGL';
% file_acquire_heat = '171212/USTplanar3M5_BA52_heat@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_29s14m00h_13-12-17_avg4_2D_raw.SGL';
%     dx = 1e-4;
%     dy = dx;
%     dt = 4e-9;
%     t_0 = 1e-6;
% file_pretune_cold = '171219/pretuning_BA52_cold/';
% file_pretune_heat = '171219/pretuning_BA52_heat/';
% file_acquire_cold  = '171219/USTplanar3M5_BA52_cold@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_12s59m15h_19-12-17_avg1_2D_raw.SGL';
% file_acquire_heat = '171219/USTplanar3M5_BA52_heat@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_32s45m15h_19-12-17_avg1_2D_raw.SGL';
%     dx = 1e-4;
%     dy = dx;
%     dt = 4e-9;
%     t_0 = 2e-6;

% BIG COMPARE EXPERIMENT: cold20Hz vs cold1kHz vs heat20Hz
% file_pretune_cold = '180105/pretuning_BA52_cold_1/';
% file_pretune_cold = '180105/pretuning_BA52_cold_2/';
% file_pretune_cold = '180105/pretuning_BA52_cold_3/';
% file_pretune_cold = '180105/pretuning_BA52_cold_4/';

% file_pretune_heat = '180105/pretuning_BA52_heat_1/';
% file_pretune_heat = '180105/pretuning_BA52_heat_2/';
% file_pretune_heat = '180105/pretuning_BA52_heat_3/';

% = '180105/USTplanar3M5_BA52_cold20Hz_1@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_28s40m11h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_cold20Hz_2@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_49s50m11h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_cold20Hz_3@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_41s56m11h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_cold20Hz_4@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_39s00m12h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_cold20Hz_5@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_07s04m12h_05-01-18_avg1_2D_raw.SGL';

% = '180105/USTplanar3M5_BA52_cold1kHz_1@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_38s07m12h_05-01-18_avg1_2D_raw.SGL';
% file_acquire_cold = '180105/USTplanar3M5_BA52_cold1kHz_2@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_55s08m12h_05-01-18_avg1_2D_raw.SGL';
% file_acquire_heat = '180105/USTplanar3M5_BA52_cold1kHz_3@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_29s09m12h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_cold1kHz_4@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_19s10m12h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_cold1kHz_5@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_41s10m12h_05-01-18_avg1_2D_raw.SGL';

% = '180105/USTplanar3M5_BA52_heat20Hz_1@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_25s06m13h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_heat20Hz_2@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_53s11m13h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_heat20Hz_3@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_15s17m13h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_heat20Hz_4@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_18s21m13h_05-01-18_avg1_2D_raw.SGL';
% = '180105/USTplanar3M5_BA52_heat20Hz_5@0nm_t0[-500]_dx[100µm]_dy[100µm]_dt[4ns]_51s25m13h_05-01-18_avg1_2D_raw.SGL';

%     dx = 1e-4;
%     dy = dx;
%     dt = 4e-9;
%     t_0 = 2e-6;

% check for centre of excitation beam
% comparePreTuningMaps(file_dir,file_pretune_cold,file_pretune_heat)

% dataSGL_cold = loadSGL( [file_dir file_acquire_cold] );
% dataSGL_heat = loadSGL( [file_dir file_acquire_heat] );

% viewing and comparing single time series
% slice_x = 30;
% slice_y = 30;
% viewSGL(dataSGL_heat(:,:,1:500),dx,dy,dt,t_0,slice_x,slice_y)

% compare sensitivity between cold and hot
% comparePeakSGL(dataSGL_cold, dataSGL_heat)

% movie comparing time series for USTplanar3M5 mappings when cold and hot
% compareSGLmovie(dataSGL_cold, dataSGL_heat, dx, dy, dt, t_0)

% compare DC during pre-tuning and acquisition
% file_params_cold = '171212/USTplanar3M5_BA52_cold@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_38s30m00h_13-12-17_avg4_2D_raw_parameters';
% file_params_heat = '171212/USTplanar3M5_BA52_heat@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_29s14m00h_13-12-17_avg4_2D_raw_parameters';

%% 180320 heatingEffects on sensitivity due to pretuning with CNT coated sensors
% find beam centre
% compare hot(20Hz), cold(20Hz), cold(1kHz)

% clear all
% close all

% file_dir = '../data/heatingEffects/180320/';

% file_dir_cold = 'findBeamCentre_pretuning_x[-3,3]_y[-3,3]_cold/';
% file_dir_heat = 'findBeamCentre_pretuning_x[-3,3]_y[-3,3]_heat/';

% file_dir_cold = 'findBeamCentre_pretuning_x[-3,3]_y[-1,5]_cold/';
% file_dir_heat = 'findBeamCentre_pretuning_x[-3,3]_y[-1,5]_heat/';

% file_dir_cold = 'findBeamCentre_pretuning_x[-5,1]_y[1,7]_cold/';
% file_dir_heat = 'findBeamCentre_pretuning_x[-5,1]_y[1,7]_heat/';

% file_dir_cold = 'findBeamCentre_pretuning_x[-1,7]_y[-1,7]_cold/';
% file_dir_heat = 'findBeamCentre_pretuning_x[-1,7]_y[-1,7]_heat/';

% beam centre
% file_dir_cold = 'findBeamCentre_pretuning_x[1,7]_y[-1,5]_cold/';
% file_dir_heat = 'findBeamCentre_pretuning_x[1,7]_y[-1,5]_heat/';

% file_dir_cold = 'PT_cold_2/';
% file_dir_heat = 'PT_heat_1/';

% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)


% file_acquire_cold = 'ACQ_BA59_cold1kHz_1@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_25s12m19h_20-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold1kHz_2@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_51s17m19h_20-03-18_avg1_2D_raw.SGL';
% bad file_acquire_cold = 'ACQ_BA59_cold1kHz_3@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_16s19m19h_20-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold1kHz_4@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_18s23m19h_20-03-18_avg1_2D_raw.SGL';

% file_acquire_cold = 'ACQ_BA59_cold20Hz_1@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_46s03m20h_20-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_2@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_31s12m20h_20-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_3@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_46s16m20h_20-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_4@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_03s21m20h_20-03-18_avg1_2D_raw.SGL';

% file_acquire_heat = 'ACQ_BA59_heat20Hz_1@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_52s36m19h_20-03-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_2@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_30s44m19h_20-03-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_3@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_57s48m19h_20-03-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_4@0nm_t0[-750]_dx[100µm]_dy[100µm]_dt[4ns]_22s53m19h_20-03-18_avg1_2D_raw.SGL';

%     dx = 1e-4;
%     dy = dx;
%     dt = 4e-9;
%     t_0 = 3e-6;

% dataSGL_cold = loadSGL( [file_dir file_acquire_cold] );
% dataSGL_heat = loadSGL( [file_dir file_acquire_heat] );

% viewing and comparing single time series
% slice_x = 30;
% slice_y = 30;
% viewSGL(dataSGL_cold(:,:,1:500),dx,dy,dt,t_0,slice_x,slice_y)
% viewSGL(dataSGL_heat(:,:,1:500),dx,dy,dt,t_0,slice_x,slice_y)

% compare sensitivity between cold and hot
% comparePeakSGL(dataSGL_cold, dataSGL_heat)


%% 180322 heatingEffects on sensitivity due to pretuning with CNT coated sensors
% find beam centre
% compare hot(20Hz), cold(20Hz), cold(1kHz)
% overlaid peak histograms from each acquisition mode
% overlaid time series from each acquisition mode

% clear all
% close all

% file_dir = '../data/heatingEffects/180322/';

% file_dir_cold = 'findBeamCentre_pretuning_x[-1,9]_y[-3,7]_cold/';
% file_dir_heat = 'findBeamCentre_pretuning_x[-1,9]_y[-3,7]_heat/';

% beam centre
% file_dir_cold = 'findBeamCentre_pretuning_x[1,7]_y[-1,5]_cold/';
% file_dir_heat = 'findBeamCentre_pretuning_x[1,7]_y[-1,5]_heat/';

% file_dir_cold = 'PT_cold_3/';
% file_dir_heat = 'PT_heat_2/';

% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)


% file_acquire_cold = 'ACQ_BA59_cold1kHz_1@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_51s08m16h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold1kHz_2@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_01s14m16h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold1kHz_3@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_30s15m16h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold1kHz_4@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_41s16m16h_22-03-18_avg1_2D_raw.SGL';

% file_acquire_cold = 'ACQ_BA59_cold20Hz_1@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_47s02m17h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_2@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_24s10m17h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_3@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_02s15m17h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_4@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_17s19m17h_22-03-18_avg1_2D_raw.SGL';

% file_acquire_heat = 'ACQ_BA59_heat20Hz_1@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_09s29m16h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_2@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_49s36m16h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_3@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_48s43m16h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_4@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_02s48m16h_22-03-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_5@0nm_t0[-2300]_dx[100µm]_dy[100µm]_dt[4ns]_44s52m16h_22-03-18_avg1_2D_raw.SGL';

%     dx = 1e-4;
%     dy = dx;
%     dt = 4e-9;
%     t_0 = 9.2e-6;

% dataSGL_cold = loadSGL( [file_dir file_acquire_cold] );
% dataSGL_heat = loadSGL( [file_dir file_acquire_heat] );

% viewing and comparing single time series
% slice_x = 30;
% slice_y = 30;
% viewSGL(dataSGL_cold(:,:,1:500),dx,dy,dt,t_0,slice_x,slice_y)
% viewSGL(dataSGL_heat(:,:,1:500),dx,dy,dt,t_0,slice_x,slice_y)

% compare sensitivity between cold and hot
% comparePeakSGL(dataSGL_cold, dataSGL_heat)

% overlaid plot histograms of peak pressures
% viewSGLpeakHisto(dataSGL_heat)
% legend('cold1kHz_1','cold1kHz_2','(cold1kHz_3)','(cold1kHz_4)')
% legend('cold20Hz_1','cold20Hz_2','cold20Hz_3','cold20Hz_4')
% legend('heat20Hz_1','heat20Hz_2','heat20Hz_3','heat20Hz_4','heat20Hz_5')


% overlaid time series from UST and laserGenUS signals
% file_name = 'waveform_laserGenUSandReflection.txt';
% file_name = 'waveform_USTplanar3M5.txt';
% file_name = 'waveform_laserGenUSandReflection_USTplanar3M5.txt';
% t_0 = -10e-6;

% file_name = 'waveformshort_laserGenUSandReflection.txt';
% file_name = 'waveformshort_USTplanar3M5[2].txt';
% file_name = 'waveformshort_laserGenUSandReflection_USTplanar3M5[2].txt';
% t_0 = 6e-6;

% viewSGLsingle(file_dir,file_name,t_0)
% legend('US transducer (cold 100Hz)','laser generated US (hot 20Hz)','UST & LGUS (hot 20 Hz)')


%% 180404 heatingEffects on sensitivity due to pretuning with CNT coated sensors
% find beam centre
% compare hot(20Hz), cold(20Hz), cold(1kHz)
% also obtain heat 20Hz with UST unplugged
% overlaid peak histograms from each acquisition mode

% clear all
% close all

% file_dir = '../data/heatingEffects/180404/';

% beam centre
% file_dir_cold = 'findBeamCentre_pretuning_x[1,7]_y[-1,5]_cold/';
% file_dir_heat = 'findBeamCentre_pretuning_x[1,7]_y[-1,5]_heat/';

% file_dir_cold = 'PT_cold_2/';
% file_dir_heat = 'PT_heat_1/';

% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)


% cold 1kHz
% file_acquire_cold = 'ACQ_BA59_cold1kHz_1@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_11s21m19h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold1kHz_2@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_45s25m19h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold1kHz_3@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_51s26m19h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold1kHz_4@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_59s27m19h_04-04-18_avg1_2D_raw.SGL';

% cold 20Hz
% file_acquire_cold = 'ACQ_BA59_cold20Hz_1@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_31s34m20h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_2@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_35s42m20h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_3@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_23s47m20h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BA59_cold20Hz_4@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_39s51m20h_04-04-18_avg1_2D_raw.SGL';

% heat 20Hz
% file_acquire_heat = 'ACQ_BA59_heat20Hz_1@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_35s39m19h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_2@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_17s54m19h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_3@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_19s03m20h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_4@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_17s12m20h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BA59_heat20Hz_5@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_40s21m20h_04-04-18_avg1_2D_raw.SGL';

% heat 20Hz with UST unplugged
% file_acquire_heat_noUST = 'ACQ_BA59_heat20Hz_1n@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_41s45m19h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_heat_noUST = 'ACQ_BA59_heat20Hz_2n@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_51s58m19h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_heat_noUST = 'ACQ_BA59_heat20Hz_3n@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_57s07m20h_04-04-18_avg1_2D_raw.SGL';
% file_acquire_heat_noUST = 'ACQ_BA59_heat20Hz_4n@0nm_t0[-2075]_dx[100µm]_dy[100µm]_dt[4ns]_38s16m20h_04-04-18_avg1_2D_raw.SGL';

%     dx = 1e-4;
%     dy = dx;
%     dt = 4e-9;
%     t_0 = 8.3e-6;

% dataSGL_cold       = loadSGL( [file_dir file_acquire_cold      ] );
% dataSGL_heat       = loadSGL( [file_dir file_acquire_heat      ] );
% dataSGL_heat_noUST = loadSGL( [file_dir file_acquire_heat_noUST] );
% dataSGL_heat_UST   = dataSGL_heat - dataSGL_heat_noUST;

% viewing and comparing single time series
% slice_x = 30;
% slice_y = 30;
% viewSGL(dataSGL_cold(:,:,1:500),dx,dy,dt,t_0,slice_x,slice_y)
% viewSGL(dataSGL_heat(:,:,1:500),dx,dy,dt,t_0,slice_x,slice_y)
% viewSGL(dataSGL_heat_UST(:,:,1:500),dx,dy,dt,t_0,slice_x,slice_y)

% compare sensitivity between cold and hot
% comparePeakSGL(dataSGL_cold, dataSGL_heat)

% overlaid plot histograms of peak pressures
% viewSGLpeakHisto(dataSGL_heat_UST)
% legend('(cold1kHz_1)','cold1kHz_2','cold1kHz_3','cold1kHz_4')
% legend('cold20Hz_1','cold20Hz_2','cold20Hz_3','cold20Hz_4')
% legend('heat20Hz_1','heat20Hz_2','heat20Hz_3','heat20Hz_4','heat20Hz_5')
% legend('heat20Hz_1 no UST','heat20Hz_2 no UST','heat20Hz_3 no UST','heat20Hz_4 no UST')
% legend('heat20Hz_1 UST','heat20Hz_2 UST','heat20Hz_3 UST','heat20Hz_4 UST')

% single time series comparison
% figure
% hold on
% plot(squeeze(dataSGL_cold(slice_x,slice_y,1:500)))
% plot(squeeze(dataSGL_heat(slice_x,slice_y,1:500)))
% plot(squeeze(dataSGL_heat_noUST(slice_x,slice_y,1:500)))
% plot(squeeze(dataSGL_heat_UST(slice_x,slice_y,1:500)))
% legend('cold20Hz','heat20Hz','heat20Hz no UST', 'heat20Hz UST')


%% 180507 heatingEffects p2p for different barrier thicknesses BA59[2.5um], BA19[9.9um], BA21[19.8um]

file_dir = '../data/heatingEffects/180507/';

file_name = 'heatingEffects_p2p_BA59[2.5um]_long_avg10_1.txt';      % 3, avg10_1
% file_name = 'heatingEffects_p2p_BA19[9.9um]_long_avg10_3.txt';      % 3, avg10_3
% file_name = 'heatingEffects_p2p_BA21[19.8um]_long_avg10_1.txt';     % 3, avg10_1
dt_sampl = 100e-9;
num_sampl = 700000;
frac_presampl = 0.1;

% file_name = 'heatingEffects_p2p_BA59[2.5um]_short_avg10_3.txt'; num_samples = 250000;
% file_name = 'heatingEffects_p2p_BA19[9.9um]_short_avg10_3.txt'; num_samples = 375000;
% file_name = 'heatingEffects_p2p_BA21[19.8um]_short_avg10_3.txt'; num_samples = 625000;
% dt_samples = 4e-9;
% frac_presamples = 0.1;

heatingEffects_p2p(file_dir,file_name,dt_sampl,num_sampl,frac_presampl,...
                'subtractBaseline',true)

% legend('BA59[2.5um] avg1','BA59[2.5um] avg10')
% legend('BA19[9.9um] avg1','BA19[9.9um] avg10')
% legend('BA21[19.8um] avg1','BA21[19.8um] avg10')
% legend('BA59[2.5um] avg10','BA19[9.9um] avg10','BA21[19.8um] avg10')
% legend('BA59[2.5um] avg10 + 0.1','BA19[9.9um] avg10','BA21[19.8um] avg10 - 0.1')

% zoom in on _long
% xlim([-0.02,0.1])
% ylim([-0.1,0.1])

% thickness = [ 2.5 9.9 19.8 ];
% delay = [ 0.04 0.3 1.1 ];
% delay_error = [ 0.01 0.1 0.1 ];
% figure
% % plot(thickness,delay,'+')
% errorbar(thickness,delay,delay_error,delay_error,'+')
%     xlabel('barrier coating thickness [um]')
%     ylabel('heat wave delay [ms]')
%     title('heat wave delay vs barrier coating thickness')
%     xlim([0,22])
%     ylim([0,1.4])


%% 180511 heatingEffects p2p for different barrier thicknesses BA59[2.5um], BA19[9.9um], BA21[19.8um]

% file_dir = '../data/heatingEffects/180511/';
% 
% sensors = {'BA59[2.5um]' 'BA19[9.9um]' 'BA21[19.8um]'}; % 
% lst_avg = [1 10 50 100];
% indeces = 1:1:3; % index: 1, 3/3, 3/2
% 
% dt_sampl = 100e-9;
% num_sampl = 200000;
% frac_presampl = 0.1;
% dt_sampl = 1e-6;
% num_sampl = 20000;
% 
% for sensor = sensors(3)
%     for num_avg = lst_avg(4)
%         for index = indeces(2)
%             file_name = ['heatingEffects_p2p_' sensor{1} '_long_avg' num2str(num_avg)...
%                             '_' num2str(index) '.txt'];
% 
%             heatingEffects_p2p(file_dir,file_name,dt_sampl,num_sampl,frac_presampl,...
%                             'subtractBaseline',true)
%         end
%     end
% end

% legend('BA59[2.5um] avg50 #1','BA59[2.5um] avg50 #2','BA59[2.5um] avg50 #3',...
%         'BA19[9.9um] avg100 #1','BA19[9.9um] avg100 #2','BA19[9.9um] avg100 #3',...
%         'BA21[19.8um] avg100 #1','BA21[19.8um] avg100 #2','BA21[19.8um] avg100 #3',...
%         'Location','northeastoutside')

% legend('BA59[2.5um] avg1','BA59[2.5um] avg10','BA59[2.5um] avg50')
% legend('BA19[9.9um] avg1','BA19[9.9um] avg10','BA19[9.9um] avg50','BA19[9.9um] avg100')
% legend('BA21[19.8um] avg1','BA21[19.8um] avg10','BA21[19.8um] avg50','BA21[19.8um] avg100')

% legend('BA59[2.5um] avg50','BA19[9.9um] avg100','BA21[19.8um] avg100')

% zoom in on _long
% xlim([-0.2,1])
% xlim([-0.5,3])
% ylim([-0.2,0.2])


%% 180813 heatingEffects on sensitivity due to pretuning with CNT coated sensors
% compare hot(20Hz), cold(20Hz), cold(1kHz)
% also obtain heat 20Hz with UST unplugged for some
% overlaid peak histograms from each acquisition mode
% 
% clear all
% close all

file_dir = '../data/heatingEffects/180813/';

% beam centre
% file_dir_cold = 'PT_cold_2/';
% file_dir_heat = 'PT_heat_1/';

% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)


% cold 1kHz
% file_acquire_cold = 'ACQ_BK31[CNT]_cold1kHz_1@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_03s56m20h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BK31[CNT]_cold1kHz_2@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_57s59m20h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BK31[CNT]_cold1kHz_3@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_54s00m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BK31[CNT]_cold1kHz_4@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_04s02m21h_13-08-18_avg1_2D_raw.SGL';

% cold 20Hz
% file_acquire_cold = 'ACQ_BK31[CNT]_cold20Hz_1@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_51s44m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BK31[CNT]_cold20Hz_2@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_05s49m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BK31[CNT]_cold20Hz_3@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_47s52m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_BK31[CNT]_cold20Hz_4@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_41s56m21h_13-08-18_avg1_2D_raw.SGL';

% heat 20Hz
% file_acquire_heat = 'ACQ_BK31[CNT]_heat20Hz_1@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_59s11m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BK31[CNT]_heat20Hz_2@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_11s21m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BK31[CNT]_heat20Hz_3@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_42s28m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BK31[CNT]_heat20Hz_4@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_39s32m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_BK31[CNT]_heat20Hz_5@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_51s37m21h_13-08-18_avg1_2D_raw.SGL';

% heat 20Hz with UST unplugged
% file_acquire_heat_noUST = 'ACQ_BK31[CNT]_heat20Hz_1n@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_56s15m21h_13-08-18_avg1_2D_raw.SGL';
% file_acquire_heat_noUST = 'ACQ_BK31[CNT]_heat20Hz_2n@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_55s24m21h_13-08-18_avg1_2D_raw.SGL';


% [dataSGL_cold , params_cold ] = loadSGL( [file_dir file_acquire_cold ] );
[dataSGL_heat , params_heat ] = loadSGL( [file_dir file_acquire_heat ] );
% dataSGL_heat_noUST = loadSGL( [file_dir file_acquire_heat_noUST] );
% dataSGL_heat_UST   = dataSGL_heat - dataSGL_heat_noUST;

% viewing and comparing single time series
t_0 = 9e-6;
slice_x = 30;
slice_y = 30;
% viewSGL(dataSGL_cold(:,:,1:400),params_cold,t_0,slice_x,slice_y)
viewSGL(dataSGL_heat(:,:,1:400),params_heat,t_0,slice_x,slice_y)
% viewSGL(dataSGL_heat_UST(:,:,1:400),dx,dy,dt,t_0,slice_x,slice_y)

% compare sensitivity between cold and hot
% comparePeakSGL(dataSGL_cold, dataSGL_heat)

% overlaid plot histograms of peak pressures
% viewSGLpeakHisto(dataSGL_heat)
% legend('cold1kHz_1','cold1kHz_2','cold1kHz_3','cold1kHz_4')
% legend('cold20Hz_1','cold20Hz_2','cold20Hz_3','cold20Hz_4')
% legend('heat20Hz_1','heat20Hz_2','heat20Hz_3','heat20Hz_4','heat20Hz_5')
% legend('heat20Hz_1 no UST','heat20Hz_2 no UST','heat20Hz_3 no UST','heat20Hz_4 no UST')
% legend('heat20Hz_1 UST','heat20Hz_2 UST','heat20Hz_3 UST','heat20Hz_4 UST')

% single time series comparison
% figure
% hold on
% plot(squeeze(dataSGL_cold(slice_x,slice_y,1:500)))
% plot(squeeze(dataSGL_heat(slice_x,slice_y,1:500)))
% plot(squeeze(dataSGL_heat_noUST(slice_x,slice_y,1:500)))
% plot(squeeze(dataSGL_heat_UST(slice_x,slice_y,1:500)))
% legend('cold20Hz','heat20Hz','heat20Hz no UST', 'heat20Hz UST')


%% checking consistency of transducer output over 2 hours

file_dir = '../data/heatingEffects/180817/';

% file_dir_cold = 'PT_1/';
% file_dir_heat = 'PT0.05_12/';
% 
% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)

% file_acquire = 'USTplanar3M5_energy2_non-inverted@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_11s54m17h_17-08-18_avg1_2D_raw.SGL';

file_acquire = 'ACQ_1@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_37s01m18h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_2@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_24s07m18h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_3@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_13s29m18h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_4@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_53s39m18h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_5@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_37s44m18h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_6@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_11s09m19h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_7@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_57s27m19h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_8@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_42s40m19h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_9@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_44s51m19h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_10@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_16s57m19h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_11@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_49s03m20h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_12@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_52s09m20h_17-08-18_avg1_2D_raw.SGL';
% file_acquire = 'ACQ_13@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_12s19m20h_17-08-18_avg1_2D_raw.SGL';

[dataSGL, params] = loadSGL([file_dir file_acquire]);
t_0 = 9.4e-6;
slice_x = 75;
slice_y = 75;

viewSGL(dataSGL(:,:,1:400),params,t_0,slice_x,slice_y)
viewSGLpeakHisto(dataSGL,'BinWidth',0.01,'Color','b')
% legend('0.1nm 1','0.1nm 2','0.1nm 3','0.1nm 4','0.1nm 5',...
%        '0.1nm 6','0.1nm 7','0.1nm 8','0.1nm 9','0.1nm 10',...
%        '0.05nm 11','0.05nm 12','0.05nm 13')


%% 180820 heatingEffects on sensitivity due to pretuning with CNT coated sensors
% compare hot(20Hz), cold(20Hz), cold(1kHz)
% overlaid peak histograms from each acquisition mode

% clear all
% close all

file_dir = '../data/heatingEffects/180820/';

% beam centre
% file_dir_cold = 'PT_cold_4/';
% file_dir_heat = 'PT_heat_4/';

% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)


% cold 1kHz
% file_acquire_cold = 'ACQ_cold1kHz_1@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_19s47m21h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_cold1kHz_2@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_55s52m21h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_cold1kHz_3@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_12s55m21h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_cold1kHz_4@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_41s57m21h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_cold1kHz_5@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_40s58m21h_20-08-18_avg1_2D_raw.SGL';

% cold 20Hz
% file_acquire_cold = 'ACQ_cold20Hz_1@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_41s40m22h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_cold20Hz_2@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_41s45m22h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_cold20Hz_3@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_32s50m22h_20-08-18_avg1_2D_raw.SGL';
file_acquire_cold = 'ACQ_cold20Hz_4@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_18s55m22h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_cold = 'ACQ_cold20Hz_5@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_52s58m22h_20-08-18_avg1_2D_raw.SGL';

% heat 20Hz
% file_acquire_heat = 'ACQ_heat20Hz_1@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_15s08m22h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_heat20Hz_2@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_17s13m22h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_heat20Hz_3@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_54s18m22h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_heat20Hz_4@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_10s24m22h_20-08-18_avg1_2D_raw.SGL';
% file_acquire_heat = 'ACQ_heat20Hz_5@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_57s27m22h_20-08-18_avg1_2D_raw.SGL';
file_acquire_heat = 'ACQ_heat20Hz_6@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[5ns]_02s32m22h_20-08-18_avg1_2D_raw.SGL';

[dataSGL_cold , params_cold ] = loadSGL( [file_dir file_acquire_cold ] );
[dataSGL_heat , params_heat ] = loadSGL( [file_dir file_acquire_heat ] );

% viewing and comparing single time series
t_0 = 9e-6;
slice_x = 30;
slice_y = 30;
% viewSGL(dataSGL_cold(:,:,1:400),params_cold,t_0,slice_x,slice_y)
% viewSGL(dataSGL_heat(:,:,1:400),params_heat,t_0,slice_x,slice_y)
% viewSGL(dataSGL_heat_UST(:,:,1:400),dx,dy,dt,t_0,slice_x,slice_y)

% compare sensitivity between cold and hot
comparePeakSGL(dataSGL_cold, dataSGL_heat)

% overlaid plot histograms of peak pressures
% viewSGLpeakHisto(dataSGL_heat)
% legend('coldPT cold1kHz_1','coldPT cold1kHz_2','coldPT cold1kHz_3','coldPT cold1kHz_4','cold1kHz_5')
% legend('coldPT cold20Hz_1','coldPT cold20Hz_2','coldPT cold20Hz_3','coldPT cold20Hz_4','cold20Hz_5')
% legend('heatPT heat20Hz_1','heatPT heat20Hz_2','heatPT heat20Hz_3','heatPT heat20Hz_4','heat20Hz_5','coldPT heat20Hz_6')


% single time series comparison
% figure
% hold on
% plot(squeeze(dataSGL_cold1kHz(slice_x,slice_y,1:400)))
% plot(squeeze(dataSGL_cold20Hz(slice_x,slice_y,1:400)))
% plot(squeeze(dataSGL_heat(slice_x,slice_y,1:400)))
% legend('cold1kHz','cold20Hz','heat20Hz')
% xlabel('time [dt]')
% ylabel('signal amplitude [V]')










