
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

