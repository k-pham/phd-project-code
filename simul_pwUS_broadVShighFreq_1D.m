% Pulse-echo plane-wave US imaging - broadband vs high frequency - 1D

clearvars;

file_dir = 'D:\PROJECT\data\simulations\freqconsi\broadVShighFreq_1D\';


%% SET UP EXPERIMENT

% create the computational grid
dx = 10e-6;                 % grid point spacing in the x direction [m]
Nx = 1024;                  % number of grid points in the x (row) direction
kgrid = kWaveGrid(Nx, dx);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define imaging object - STEP
% step_x = round(Nx/2);
% step_c = 1600;
% step_rho = 1100;

% define imaging object - RAMP
ramp_width = 10;
ramp_xmin = round(Nx/2) - ramp_width/2;
ramp_xmax = round(Nx/2) + ramp_width/2 - 1;
ramp_x = ramp_xmin : ramp_xmax;
ramp_cmax = 1600;
ramp_rhomax = 1100;
ramp_c = linspace(c0, ramp_cmax, ramp_width);
ramp_rho = linspace(rho0, ramp_rhomax, ramp_width);

% define sound speed and density of medium
medium.sound_speed = c0 * ones(Nx, 1);
medium.density = rho0 * ones(Nx, 1);

% medium.sound_speed(step_x : end) = step_c;
% medium.density(step_x : end) = step_rho;
medium.sound_speed(ramp_x) = ramp_c;
medium.sound_speed(ramp_xmax+1:end) = ramp_cmax;
medium.density(ramp_x) = ramp_rho;
medium.density(ramp_xmax+1:end) = ramp_rhomax;

% create time array
cfl = 0.2;
t_end = 2*Nx*dx/c0;
kgrid.makeTime(medium.sound_speed,cfl,t_end);

% define source
source.p_mask = zeros(Nx, 1);
source.p_mask(pml_size + 1) = 1;
source.p = gaussian(kgrid.t_array,1,20e-9,(6.8e-9)^2);

% define sensor position at source
sensor.mask = source.p_mask;


%% SIMULATE EXPERIMENT

inputs = {'PMLSize', pml_size, 'PlotLayout', true, 'PlotSim', true};
sensor_data = kspaceFirstOrder1D(kgrid, medium, source, sensor, inputs{:});

save([file_dir 'sensor_data_1D'], 'sensor_data')

figure
plot(kgrid.t_array,sensor_data)
ylabel('amplitude')
xlabel('time / \mus')


















