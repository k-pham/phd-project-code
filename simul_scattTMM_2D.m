% Pulse-echo plane-wave US imaging - scattering TMM with holes (simulate agarTMM)
% combine _rand and _points scripts into one

global kgrid t_array

%% loops for parameter search

for c_range = 0:10:150
for rho_range = 0:10:100
%% SET UP

% clear all
close all

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';


%% SET UP EXPERIMENT

% create the computational grid
dx = 10e-6;                 % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
Nx = 1536;                  % number of grid points in the x (row) direction
Ny = 1024;                  % number of grid points in the y (column) direction
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the thickness of the PML [grid points]
pml_size = 20;


%% MEDIUM

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% choose option to represent scattering medium as
scattering_type = 'random';      % options: 'random', 'points'

% define scattering medium as random medium
% c_range   = 0;
% rho_range = 0;

% define scattering medium as point scatterers
c_pointscatt   = 1550;
rho_pointscatt = 1100;
    num_points_per_voxel = 2;
    vox_size = 100e-6;

% define non-scattering holes
c_hole      = 1500;
rho_hole    = 1000;

switch scattering_type
    case 'random'
        medium = define_random_medium(Nx,Ny,c0,rho0,c_range,rho_range,c_hole,rho_hole);
        c_scatt = c_range;
        rho_scatt = rho_range;
    case 'points'
        medium = define_pointscatt_medium(Nx,Ny,c0,rho0,c_pointscatt,rho_pointscatt,dx,dy,num_points_per_voxel,vox_size,c_hole,rho_hole);
        c_scatt = c_pointscatt;
        rho_scatt = rho_pointscatt;
end

% plot medium sound speed and density
fig_medium = figure;
set(gcf,'Position',[200,200,500,700])
subplot(2,1,1)
imagesc(kgrid.x_vec*1e3,kgrid.y_vec*1e3,medium.sound_speed')
    axis image
    title(['sound speed: ' scattering_type num2str(c_scatt)])
    xlabel('x position / mm')
    ylabel('y position / mm')
    colorbar
subplot(2,1,2)
imagesc(kgrid.x_vec*1e3,kgrid.y_vec*1e3,medium.density')
    axis image
    title(['density: ' scattering_type num2str(rho_scatt)])
    xlabel('x position / mm')
    ylabel('y position / mm')
    colorbar


%% TIME

cfl = 0.2;              % CFL number
% t_end = 2*Ny*dy/c0;     % end time of the simulation [s]
t_end = 1.5*Ny*dy/c0;     % end time of the simulation [s]
kgrid.makeTime(medium.sound_speed,cfl,t_end);


%% SOURCE

source_width = 0.4;         % proportion of width Nx
source_amplitude = 10;      % [Pa}

source.p_mask = zeros(Nx, Ny);
source.p_mask(:, pml_size + 1) = 1;
pressure = gaussian(kgrid.t_array,1,20e-9,(6.8e-9)^2);
apodisation = getWin(Nx, 'Gaussian', 'Param', source_width);
source.p = source_amplitude * apodisation * pressure;


%% SENSOR

% sensor array co-aligned with the source
sensor_positions_x = pml_size:10:(Nx - pml_size);
sensor.mask = zeros(Nx, Ny);
sensor.mask(sensor_positions_x, pml_size+1) = 1;

% calculate the number of sensor elements used
num_sensors = sum(sensor.mask(:));


%% SIMULATE EXPERIMENT

inputs = {'PMLSize', pml_size, 'PlotLayout', true, 'PlotSim', true};
sensor_data = kspaceFirstOrder2DC(kgrid, medium, source, sensor, inputs{:});

% save simulation set-up
simu.kgrid  = kgrid;
simu.medium = medium;
simu.source = source;
simu.sensor = sensor;
simu.inputs = inputs;

% figure
% imagesc(sensor_data(:,1:50)')
%     xlabel('x position / dx')
%     ylabel('time / dt')

fig_data = figure;
imagesc(kgrid.x_vec*1e3,kgrid.t_array(50:end)*1e6,sensor_data(:,50:end)')
    title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt)])
    xlabel('x position / mm')
    ylabel('time / \mus')
    colorbar


%% IMAGE FORMATION

params.Nx = size(sensor_data,1);
params.Ny = 1;
params.dx = 10*dx;
params.dy = dy;
params.dt = kgrid.dt;

params.trigger_delay        = 0;
params.Nt_zero_pad_source   = 50;
params.Nt_t0_correct        = -16;
params.file_data            = '111111\scattTMM_simul';

reflection_image = reconstruct2dUSimage(sensor_data, params, c0);
    % NOTE: kgrid and t_array UPDATED

fig_image = figure;
% imagesc(kgrid.x_vec*1e3,kgrid.t_array*c0/2*1e3,reflection_image')
%     % can use kgrid.x_vec here, even though spacing is different in image,
%     % because only care about end points of sensor; kgrid.t_array is
%     % correct to use with c0/2 scaling
imagesc(kgrid.x_vec*1e3,t_array*c0*1e3,reflection_image')
    % using updated kgrid and t_array from reconstruct2dUSimage
    axis image
    title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt)])
    xlabel('x position / mm')
    ylabel('y position / mm')
    colorbar

% prepare for saving
[Nx_image, Ny_image] = size(reflection_image);
volume_data = reshape(reflection_image,Nx_image,1,Ny_image);
volume_spacing = [kgrid.dx, 1, params.dt*c0];           % omit factor 1/2 in dz because of doubled depth bug


%% SAVE FIGURES & SENSOR DATA & IMAGE DATA

file_name = [scattering_type '_SCATT_c' num2str(c_scatt) '_rho' num2str(rho_scatt) ...
                             '_HOLE_c' num2str(c_hole) '_rho' num2str(rho_hole) ];

% saveas(fig_medium,[file_dir_figs file_name '_medium.fig'])
saveas(fig_medium,[file_dir_figs file_name '_medium.jpg'])
% saveas(fig_data,  [file_dir_figs file_name '_data.fig'])
saveas(fig_data,  [file_dir_figs file_name '_data.jpg'])
% saveas(fig_image, [file_dir_figs file_name '_image.fig'])
saveas(fig_image, [file_dir_figs file_name '_image.jpg'])

save([file_dir_data file_name '_sensor_data'], 'sensor_data', 'params', 'simu')

save([file_dir_data file_name '_image_data.mat'], 'volume_data', 'volume_spacing', 'kgrid', 't_array', '-v7.3')


%% end of loops for parameter search
end
end

%% FUNCTIONS

function medium = define_random_medium(Nx,Ny,c0,rho0,c_range,rho_range,c_hole,rho_hole)

    c_rand    = (rand(Nx,Ny)-0.5) * c_range;
    rho_rand  = (rand(Nx,Ny)-0.5) * rho_range;
    holes = get_hole_locations(Nx, Ny);
    
    % define sound speed and density of medium
    medium.sound_speed = c0   * ones(Nx, Ny) + c_rand;
    medium.density     = rho0 * ones(Nx, Ny) + rho_rand;

    medium.sound_speed(holes==1) = c_hole;
    medium.density(holes==1)     = rho_hole;

end

function medium = define_pointscatt_medium(Nx,Ny,c0,rho0,c_pointscatt,rho_pointscatt,dx,dy,num_points_per_voxel,vox_size,c_hole,rho_hole)

    pointscatts = get_pointscatt_locations(Nx, Ny, dx, dy, num_points_per_voxel, vox_size);
    holes = get_hole_locations(Nx, Ny);
    
    % define sound speed and density of medium
    medium.sound_speed = c0   * ones(Nx, Ny);
    medium.density     = rho0 * ones(Nx, Ny);

    medium.sound_speed(pointscatts==1) = c_pointscatt;
    medium.density(pointscatts==1)     = rho_pointscatt;

    medium.sound_speed(holes==1) = c_hole;
    medium.density(holes==1)     = rho_hole;

end

function pointscatts = get_pointscatt_locations(Nx, Ny, dx, dy, num_points_per_voxel, vox_size)

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

end

function holes = get_hole_locations(Nx, Ny)

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
    
end

function slab = get_slab_location(Nx, Ny)

    slab_thickness = 100;
    slab_y         = 250;
    slab_yrange    = round(slab_y-slab_thickness/2+0.5 : slab_y+slab_thickness/2-0.5);
    
    slab = zeros(Nx,Ny);
    slab(:,slab_yrange) = 1;
    
end
