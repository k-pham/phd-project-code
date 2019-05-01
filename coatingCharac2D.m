clear all
close all

%% Carbon nanotubes on BA58 (paths)

file_dir = '../data/coatingCharac/';
% file_name_p2p = '171101/heatingEffects_BA58[nanoparticle]_pulse2pulse.txt';
%     dt_samples_p2p = 100e-9;
%     num_samples_p2p = 700000;
% 
% file_dir_cold = '171101/heatingEffects_BA58[nanoparticle]_pretuning_cold/';
% file_dir_heat = '171101/heatingEffects_BA58[nanoparticle]_pretuning_heating/';

%file_name_laserGenUS = '171107/ULTRA3_143us_nanoparticle_harddielectric@0nm_t0_-2250__dx_100µm__dy_100µm__dt_4ns__48s50m19h_07-11-17_avg1_2D_raw.SGL';
%file_name_laserGenUS = '171114/ULTRA3[143us]_diffuser_nanoparticle_harddielectric@0nm_t0[-2250]_dx[100µm]_dy[100µm]_dt[4ns]_27s09m17h_14-11-17_avg1_2D_raw.SGL';
%file_name_laserGenUS = '171114/ULTRA3[143us]_diffuser_nanoparticle_harddielectric@0nm_t0[-2250]_dx[100µm]_dy[100µm]_dt[4ns]_27s09m17h_14-11-17_avg1_2D_corr.SGL';
%file_name_laserGenUS = '171114/ULTRA3[143us]_diffuser_spraypaint_harddielectric@0nm_t0[-2250]_dx[100µm]_dy[100µm]_dt[4ns]_14s22m17h_14-11-17_avg1_2D_raw.SGL';
%file_name_laserGenUS = '171117/ULTRA3[143us]_diffuser_nanoparticle_harddielectric@0nm_t0[-2250]_dx[100µm]_dy[100µm]_dt[4ns]_19s58m19h_17-11-17_avg4_2D_raw.SGL'; %moved excitation after
%file_name_laserGenUS = '171117/ULTRA3[143us]_diffuser_spraypaint_harddielectric@0nm_t0[-2250]_dx[100µm]_dy[100µm]_dt[4ns]_17s39m20h_17-11-17_avg4_2D_raw.SGL';
%file_name_laserGenUS = '171117/ULTRA3[143us]_diffuser_nanoparticle_harddielectric@0nm_t0[-2250]_dx[100µm]_dy[100µm]_dt[4ns]_34s22m21h_17-11-17_avg4_2D_raw.SGL';
%file_name_laserGenUS = '171205/ULTRA3[143us]_diffuser2_spraypaint_harddielectric@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_17s09m22h_05-12-17_avg8_2D_raw.SGL';
%file_name_laserGenUS = '171205/ULTRA3[143us]_diffuser2_nanotube_harddielectric@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_54s00m00h_06-12-17_avg8_2D_raw.SGL';
%file_name_laserGenUS = '171211/ULTRA3[143us]_diffuser2ii_ink[acrylic]_harddielectric@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_05s13m21h_11-12-17_avg8_2D_raw.SGL';
%file_name_laserGenUS = '171212/ULTRA3[143us]_diffuser2ii_persplex_ink[acrylic]_harddielectric@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_29s12m21h_12-12-17_avg12_2D_raw.SGL';
%file_name_laserGenUS = '171212/ULTRA3[143us]_diffuser2ii_nanotubes_harddielectric@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_15s29m22h_12-12-17_avg4_2D_raw.SGL';

%file_name_laserGenUS = '180116/ULTRA3[143us]_diffuser1b_CNT_allharddielectric@0nm_t0[-875]_dx[100µm]_dy[100µm]_dt[4ns]_33s51m16h_16-01-18_avg1_2D_raw.SGL';
%file_name_laserGenUS = '180116/ULTRA3[143us]_diffuser1b_SP_allharddielectric@0nm_t0[-875]_dx[100µm]_dy[100µm]_dt[4ns]_29s23m17h_16-01-18_avg1_2D_raw.SGL';
% file_name_laserGenUS = '180201/ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric@0nm_t0[-1250]_dx[100µm]_dy[100µm]_dt[4ns]_51s01m16h_01-02-18_avg1_2D_raw.SGL';
%file_name_laserGenUS = '180201/ULTRA3[143us]_diffuser1b05_SP[perspex]_allharddielectric@0nm_t0[-1250]_dx[100µm]_dy[100µm]_dt[4ns]_27s20m16h_01-02-18_avg1_2D_raw.SGL';
%file_name_laserGenUS = '180206/ULTRA3[143us]_diffuser1b05_CNT[BA59]_allharddielectric@0nm_t0[-985]_dx[100µm]_dy[100µm]_dt[4ns]_53s08m22h_06-02-18_avg1_2D_raw.SGL';
%     dx = 1e-4;
%     dy = dx;
%     dt = 4e-9;
%     t_0 = 5e-6;
    
    
% file_name = '181003/OPO600_AuNP[perspex]_AHD1@0nm_t0[-625]_dx[100µm]_dy[100µm]_dt[4ns]_04s37m14h_03-10-18_avg1_2D_raw.SGL';
%     t_0 = 2e-6;
    
% file_name = '181012/AuNPfilm_AHD1@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_49s21m18h_12-10-18_avg1_2D_raw.SGL';
% file_name = '181012/AuNPfilm_flipped_AHD1@0nm_t0[60]_dx[100µm]_dy[100µm]_dt[4ns]_57s10m18h_12-10-18_avg1_2D_raw.SGL';
% file_name = '181023/OPO600_AuNP[unbacked]_AHD1_avg100@0nm_t0[0]_dx[0µm]_dy[100µm]_dt[4ns]_29s35m12h_23-10-18_avg100_1D_raw.SGL';
% file_name = '181023/OPO600_AuNP[unbacked]_flipped_AHD1_avg100@0nm_t0[0]_dx[0µm]_dy[100µm]_dt[4ns]_49s50m12h_23-10-18_avg100_1D_raw.SGL';
%     t_0 = 0;

file_name = '190411/ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_wholefield@0nm_t0[0]_dx[200µm]_dy[200µm]_dt[4ns]_16s59m18h_11-04-19_avg1_2D_raw.SGL';
    t_0 = 0;
file_name = '190411/ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_wavefront@0nm_t0[-250]_dx[100µm]_dy[100µm]_dt[4ns]_27s38m19h_11-04-19_avg1_2D_raw.SGL';
    t_0 = 2e-6;
file_name = '190411/ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_31s46m19h_11-04-19_avg1_2D_raw.SGL';
    t_0 = 0;

%% heatingEffects: pulse2pulse
% heatingEffects_p2p(file_dir,file_name_p2p,dt_samples_p2p,num_samples_p2p)

%% heatingEffects: pre-tuning maps cold/hot comparison
% comparePreTuningMaps(file_dir,file_dir_cold,file_dir_heat)

%% viewSGL: map laserGenUS in transmission  &  freqSpecSGL
[dataSGL, params] = loadSGL( [file_dir file_name] );
slice_x = 4;
slice_y = 4;
viewSGL(dataSGL,params,t_0,slice_x,slice_y)
% for slice = 1:150
%     viewSGLline(dataSGL,params,t_0,slice)
% %     pause
% end
% t_min = 100;
% t_max = 300;
% freqSpecSGL(dataSGL,1/dt,slice_x,slice_y,t_min,t_max)
% [frequency, f_series_avg] = freqSpecSGLavgLine(dataSGL(:,700:999),params);

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

%% uniformity
%histogram plot of PMAX
figure
histogram(PMAX_vec(mask),'BinWidth',0.02)
    title('histogram of peak pressure across scanning area')
    xlabel('peak signal amplitude')
    ylabel(['count (out of ' num2str(Nx*Ny) ')'])
    %legend('whole scanning area','within ROI')

%% non-linear propagation?
%scatter plotting PMAX vs TOA   -   NOT APPROPRIATE
figure
scatter(TOA_vec(mask),PMAX_vec(mask),'.')
hold on
    title('correlation between peak pressure amplitude and TOA')
    xlabel('TOA index')
    ylabel('peak signal amplitude')
    axis([150,200,-0.1,0.7])
    %legend('whole scanning area','within ROI')

%% non-linear propagation? v2.0
%scatter plotting PMAX vs pulse half width ratio   -   NOT APPROPRIATE

%threshold dataSGL at half PMAX for every interrogation point
dataSGLthres_xy = zeros(size(dataSGL));
TOA_rise_xy = zeros([Nx Ny]);
TOA_fall_xy = zeros([Nx Ny]);
for xindex = 1:Nx
   for yindex = 1:Ny
      PMAX = PMAX_xy(xindex,yindex);
      if PMAX < 0.32
          continue
      end
      dataSGLthres = squeeze(-dataSGL(xindex,yindex,:));
      dataSGLthres(dataSGLthres > PMAX/2) = 1;
      dataSGLthres(dataSGLthres < PMAX/2) = 0;
      dataSGLthres_xy(xindex,yindex,:) = dataSGLthres;
      
      %find TOA of rising/falling edge of pulses from thresholded dataSGL
      TOA_prerise = find(dataSGLthres,1) - 1;
      TOA_prefall = find(dataSGLthres,1,'last');
      
      TOA_rise = interp1( [-dataSGL(xindex,yindex,TOA_prerise) -dataSGL(xindex,yindex,TOA_prerise+1)] , ...
                          [TOA_prerise TOA_prerise+1 ], ...
                           PMAX/2, ...
                          'linear' );
	  
      TOA_fall = interp1( [-dataSGL(xindex,yindex,TOA_prefall) -dataSGL(xindex,yindex,TOA_prefall+1)] , ...
                          [TOA_prefall TOA_prefall+1 ], ...
                           PMAX/2, ...
                          'linear' );
      
      TOA_rise_xy(xindex,yindex) = TOA_rise;
      TOA_fall_xy(xindex,yindex) = TOA_fall;
      
   end
end

pulse_shape_xy = (TOA_xy - TOA_rise_xy) ./ (TOA_fall_xy - TOA_xy);
pulse_shape_vec  = reshape(pulse_shape_xy, [Nx*Ny 1]);

xindex = 80;
yindex = 62;
figure('pos',[250 150 500 200])
plot(squeeze(-dataSGL(xindex,yindex,:)))
hold on
plot(squeeze(dataSGLthres_xy(xindex,yindex,:)))
plot(TOA_rise_xy(xindex,yindex),PMAX_xy(xindex,yindex)/2,'rx')
plot(TOA_xy(xindex,yindex),PMAX_xy(xindex,yindex),'rx')
plot(TOA_fall_xy(xindex,yindex),PMAX_xy(xindex,yindex)/2,'rx')
    axis([162,182,-0.1,1.1])
    title(['time series zoomed in on pulse and thresholding (x = ',num2str(xindex),', y = ',num2str(yindex), ')'])
    xlabel('TOA index')
    ylabel('signal amplitude')

figure
scatter(pulse_shape_vec,PMAX_vec,'.')
hold on
    title('correlation between peak pressure amplitude and pulse shape (half width ratio)')
    xlabel('pulse shape (half width ratio)')
    ylabel('peak signal amplitude')
    axis([0,1.2,0.2,1.4])
    %axis([0,1.2,0.2,0.8])
    %legend('whole scanning area','within ROI')

%% planarity plots
%3d scatter plot TOA over 2d grid
figure
scatter3(Xv,Yv,TOA_vec,'.')
    title('Time index of arrival for each interrogation point')
    xlabel('x')
    ylabel('y')
    zlabel('time index of arrival')
    zlim([165,185])

%2d plot of TOA in colour
figure
imagesc(TOA_xy)
    title('Time index of arrival for each interrogation point')
    xlabel('y')
    ylabel('x')
    zlabel('time index of arrival')
    colorbar
    %caxis([224 234])

%% planarity fit
%get mask for planarity
% mask = roipoly;

%apply mask to data to use for regression
TOA_vec_mask = TOA_xy(mask);
Xv_mask = Xv(mask);
Yv_mask = Yv(mask);

%regression in 2d TOA(x,y) = p1 + p2*x + p3*y using masked data points
Model = [ones(length(Xv_mask),1) Xv_mask Yv_mask];
params = Model \ TOA_vec_mask;

%evaluate fit over whole scanning area
TOA_fit_vec = params(1) + params(2) * Xv + params(3) * Yv;
TOA_fit_xy = reshape(TOA_fit_vec,[Nx Ny]);

%get deviation in time index of arrival from fit within mask
TOA_fit_vec_mask = TOA_fit_vec(mask);
deviation_mask = TOA_vec_mask - TOA_fit_vec_mask;

%plot deviation of TOA from TOA_fit in 2d colour scatter
figure
scatter(Xv_mask,Yv_mask,[],deviation_mask,'.')
hold on
    colorbar
    caxis([-3,3])
    title('time index of arrival: deviation from fit within ROI')
    xlabel('x')
    ylabel('y')
    xlim([30,120])
    ylim([20,120])
    axis equal
    
%plot 3d scatter of TOA for whole scanning area with fit from masked data
figure
scatter3(Xv,Yv,TOA_vec,'.')
hold on
mesh(X,Y,TOA_fit_xy)
    title('Time index of arrival for each interrogation point with fit')
    xlabel('x')
    ylabel('y')
    zlabel('time index of arrival')
    xlim([0,Nx])
    ylim([0,Ny])
    zlim([160,180])

