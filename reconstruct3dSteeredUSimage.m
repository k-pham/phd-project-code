%% load data

file_dir  = 'D:\PROJECT\data\angleComp\210630\';
file_data = 'polymerLeaf_angledCNT_0@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_26s35m20h_30-06-21_avg1_2D_raw.SGL';

[sensor_data, params] = loadSGL([file_dir file_data]);

fig_data = figure;
set(gcf,'Position',[100,50,600,800])
imagesc(squeeze(sensor_data(70,:,:))')
    colormap(gray)
    colorbar
    title(strtok(file_data,'@'),'Interpreter','None')
    xlabel('x axis [dx]')
    ylabel('time [dt]')
    drawnow
fig_data2 = figure;
set(gcf,'Position',[500,50,800,600])
plot(squeeze(sensor_data(70,70,:)))
    xlabel('time [dt]')
    ylabel('x axis [dx]')


%% fit plane wave to source

% note peaks will be saturated, so for TOA should pick first one at max val

% get PMAX and TOA in 2d grid and vector
[ PMAX_xy, TOA_xy ] = max(sensor_data,[],3);  % if max element occurs more than once, then TOA contains index to first occurrence
[Nx,Ny] = size(TOA_xy);
x = 1:Nx;
y = 1:Ny;
[X,Y] = meshgrid(x,y);X=X';Y=Y';
Xv = reshape(X,[Nx*Ny 1]);
Yv = reshape(Y,[Nx*Ny 1]);
PMAX_vec = reshape(PMAX_xy,[Nx*Ny 1]);
TOA_vec  = reshape(TOA_xy, [Nx*Ny 1]);

% may need to apply mask for central area to avoid outliers at edge ruining fit

% % apply mask to data to use for regression
% TOA_vec_mask = TOA_xy(mask);
% Xv_mask = Xv(mask);
% Yv_mask = Yv(mask);
TOA_vec_mask = TOA_vec;
Xv_mask = Xv;
Yv_mask = Yv;

% regression in 2d TOA(x,y) = p1 + p2*x + p3*y using masked data points
Model = [ones(length(Xv_mask),1) Xv_mask Yv_mask];
params = Model \ TOA_vec_mask;

% evaluate fit over whole scanning area
TOA_fit_vec = params(1) + params(2) * Xv + params(3) * Yv;
TOA_fit_xy = reshape(TOA_fit_vec,[Nx Ny]);

% get deviation in time index of arrival from fit within mask
TOA_fit_vec_mask = TOA_fit_vec;%(mask);
deviation_mask = TOA_vec_mask - TOA_fit_vec_mask;

% plot deviation of TOA from TOA_fit in 2d colour scatter
figure
scatter(Xv_mask,Yv_mask,[],deviation_mask,'o')
hold on
    colorbar
    %caxis([-3,3])
    title('time index of arrival: deviation from fit within ROI')
    xlabel('x')
    ylabel('y')
    %xlim([30,120])
    %ylim([20,120])
    axis equal
    
% plot 3d scatter of TOA for whole scanning area with fit from masked data
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
    %zlim([160,180])


%% calculate steering angle (a), rotation angle (b) & time offset (t0)



%% t0 correction



%% zero padding delay


