% Pulse-echo plane-wave US imaging - broadband vs high frequency

clearvars;

file_dir = 'D:\PROJECT\data\simulations\freqconsi\slowImpChange_2D\';


%% SET UP EXPERIMENT

% create the computational grid
dx = 10e-6;                 % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
Nx = 1024;                  % number of grid points in the x (row) direction
Ny = 1024;                  % number of grid points in the y (column) direction
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define imaging target
% blob_radius = 50;
% blob_x = round(Nx/2);
% blob_y = round(Ny/2);
% blob = makeDisc(Nx, Ny, blob_x, blob_y, blob_radius);
% blob_c = repmat(linspace(1500,1700,Ny),Nx,1);
% blob_rho = repmat(linspace(1200,1600,Ny),Nx,1);



% define sound speed and density of medium
medium.sound_speed = c0 * ones(Nx, Ny);          % sound speed [m/s]
medium.density = rho0 * ones(Nx, Ny);            % density [kg/m^3]

medium.sound_speed(blob == 1) = blob_c(blob == 1);
medium.density(blob == 1) = blob_rho(blob == 1);

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
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputs{:});

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
params.Nt_zero_pad_source   = 500;
params.Nt_t0_correct        = -262;
params.file_data            = '111111\broadVShighFreq_simul';

reflection_image = reconstruct2dUSimage(sensor_data, params, c0);






