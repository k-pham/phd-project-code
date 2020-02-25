% Pulse-echo plane-wave US imaging - scattering TMM with holes (simulate agarTMM)
% simulate scattering as randomly heterogeneous medium

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

% define scattering medium and non-scattering holes
c_range   = 0;
c_rand    = (rand(Nx,Ny)-0.5) * c_range;
rho_range = 100;
rho_rand  = (rand(Nx,Ny)-0.5) * rho_range;

num_holes   = 4;
hole_radius = 50;
hole_xs     = round((1:1:num_holes)*Nx/(num_holes+1));
hole_ys     = round((1:1:num_holes)*Ny/(num_holes+1));
holes = zeros(Nx,Ny);
for i = 1:num_holes
    holes = holes + makeDisc(Nx,Ny,hole_xs(i),hole_ys(i),hole_radius);
end

% define sound speed and density of medium
medium.sound_speed = c0   * ones(Nx, Ny) + c_rand;
medium.density     = rho0 * ones(Nx, Ny) + rho_rand;

medium.sound_speed(holes==1) = c0;
medium.density(holes==1) = rho0;

% create time array
cfl = 0.2;              % CFL number
t_end = 2*Ny*dy/c0;     % end time of the simulation [s]
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
imagesc(kgrid.x_vec*1e3,kgrid.t_array*1e6,sensor_data')
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

figure
imagesc(kgrid.x_vec*1e3,kgrid.y_vec*1e3,reflection_image')
xlabel('x position / mm')
ylabel('y position / mm')




