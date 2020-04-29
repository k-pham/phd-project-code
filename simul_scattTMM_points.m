% Pulse-echo plane-wave US imaging - scattering TMM with holes (simulate agarTMM)
% simulate scattering as point scatterers in otherwise homogeneous medium

clear all
close all

file_dir = 'D:\PROJECT\data\simulations\scattTMM\';


%% SET UP EXPERIMENT

% create the computational grid
dx = 10e-6;                 % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
Nx = 1536;                  % number of grid points in the x (row) direction
Ny = 1024;                  % number of grid points in the y (column) direction
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define scattering medium
c_pointscatt   = 1550;
rho_pointscatt = 1100;

num_points_per_voxel = 2;
vox_size = 100e-6;
num_vox_x = floor(Nx*dx/vox_size);
num_vox_y = floor(Ny*dy/vox_size);
vox_x_bounds = round(linspace(1,Nx+1,num_vox_x+1));
vox_y_bounds = round(linspace(1,Ny+1,num_vox_y+1));

pointscatts = zeros(Nx,Ny);
for idx_x = 1 : num_vox_x
    x_left  = vox_x_bounds(idx_x);
    x_right = vox_x_bounds(idx_x+1)-1;
    for idx_y = 1 : num_vox_y
        y_top    = vox_y_bounds(idx_y);
        y_bottom = vox_y_bounds(idx_y+1)-1;
        
        current_vox_pointscatts = zeros(x_right-x_left+1, y_bottom-y_top+1);
        scatt_locations = randsample(numel(current_vox_pointscatts), num_points_per_voxel);
        current_vox_pointscatts(scatt_locations) = 1;
        
        pointscatts(x_left:x_right,y_top:y_bottom) = current_vox_pointscatts;
        
    end
end

% define non-scattering holes
c_hole      = 1500;
rho_hole    = 1300;

% num_holes   = 4;
num_holes   = 1;
hole_radius = 50;
% hole_xs     = round((1:1:num_holes)*Nx/(num_holes+1));
% hole_ys     = round((1:1:num_holes)*Ny/(num_holes+1));
hole_xs     = round(Nx/2);
hole_ys     = round(Ny/4);
holes = zeros(Nx,Ny);
for i = 1:num_holes
    holes = holes + makeDisc(Nx,Ny,hole_xs(i),hole_ys(i),hole_radius);
end

% define sound speed and density of medium
medium.sound_speed = c0   * ones(Nx, Ny);
medium.density     = rho0 * ones(Nx, Ny);

medium.sound_speed(pointscatts==1) = c_pointscatt;
medium.density(pointscatts==1)     = rho_pointscatt;

medium.sound_speed(holes==1) = c_hole;
medium.density(holes==1)     = rho_hole;

% plot medium sound speed and density
figure
set(gcf,'Position',[200,200,500,700])
subplot(2,1,1)
imagesc(kgrid.x_vec*1e3,kgrid.y_vec*1e3,medium.sound_speed')
    axis image
    title('sound speed')
    xlabel('x position / mm')
    ylabel('y position / mm')
    colorbar
subplot(2,1,2)
imagesc(kgrid.x_vec*1e3,kgrid.y_vec*1e3,medium.density')
    axis image
    title('density')
    xlabel('x position / mm')
    ylabel('y position / mm')
    colorbar

% create time array
cfl = 0.2;              % CFL number
% t_end = 2*Ny*dy/c0;     % end time of the simulation [s]
t_end = 1.5*Ny*dy/c0;     % end time of the simulation [s]
kgrid.makeTime(medium.sound_speed,cfl,t_end);

% define source
source.p_mask = zeros(Nx, Ny);
source.p_mask(:, pml_size + 1) = 1;
sourcewidth = 0.4;
pressure = gaussian(kgrid.t_array,1,20e-9,(6.8e-9)^2);
apodisation = getWin(Nx, 'Gaussian', 'Param', sourcewidth);
source.p = apodisation * pressure;

% define a sensor array co-aligned with the source
sensor_positions_x = pml_size:10:(Nx - pml_size);
sensor.mask = zeros(Nx, Ny);
sensor.mask(sensor_positions_x, pml_size+1) = 1;

% calculate the number of sensor elements used
num_sensors = sum(sensor.mask(:));


%% SIMULATE EXPERIMENT

inputs = {'PMLSize', pml_size, 'PlotLayout', true, 'PlotSim', true};
sensor_data = kspaceFirstOrder2DC(kgrid, medium, source, sensor, inputs{:});

save([file_dir 'sensor_data'], 'sensor_data')

figure
imagesc(kgrid.x_vec*1e3,kgrid.t_array(50:end)*1e6,sensor_data(:,50:end)')
    xlabel('x position / mm')
    ylabel('time / \mus')


%% IMAGE FORMATION

params.Nx = size(sensor_data,1);
params.Ny = 1;
params.dx = 10*dx;
params.dy = dy;
params.dt = kgrid.dt;

params.trigger_delay        = 0;
params.Nt_zero_pad_source   = 50;
params.Nt_t0_correct        = -17;
params.file_data            = '111111\scattTMM_simul';

reflection_image = reconstruct2dUSimage(sensor_data, params, c0);

fig_img = figure;
imagesc(kgrid.x_vec*1e3,kgrid.y_vec*1e3,reflection_image')
    axis image
    xlabel('x position / mm')
    ylabel('y position / mm')


%% IMAGE QUALITY ASSESSMENT

% % outline of holes
% holes_outline = zeros(Nx,Ny);
% for i = 1:num_holes
%     holes_outline = holes_outline + makeCircle(Nx,Ny,hole_xs(i),hole_ys(i),hole_radius);
% end
% figure
% % hold on
% imagesc(kgrid.x_vec*1e3,kgrid.y_vec*1e3,holes_outline')

signal_hole = get_peak_signal_of_hole(reflection_image);
[scatter_hole_mean, scatter_hole_std] = get_scattering_distr_in_hole(reflection_image);



function signal_hole = get_peak_signal_of_hole(reflection_image) % hard numbers !
% limits customised to third tube

    signal_hole_front = max(max(reflection_image(116:125,1700:1900)));
    signal_hole_back  = max(max(reflection_image(116:125,2250:2450)));
    signal_hole = mean([signal_hole_front,signal_hole_back]);

end

function [scatter_hole_mean, scatter_hole_std] = get_scattering_distr_in_hole(reflection_image) % hard numbers

    ROI = reflection_image(116:125,1950:2200);
    
    scatter_hole_mean = mean(ROI(:));
    scatter_hole_std  = std(ROI(:));

end
