% Pulse-echo plane-wave US imaging - broadband vs high frequency

clearvars;

file_dir = 'D:\PROJECT\data\simulations\freqbandwidth\';


%% SET UP EXPERIMENT

% create the computational grid
dx = 10e-6;                 % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
Nx = 1024;                  % number of grid points in the x (row) direction
Ny = 700;                  % number of grid points in the y (column) direction
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define imaging target
blob_radius = 50;
blob_x = 512;
blob_y = 512;
blob = makeDisc(Nx, Ny, blob_x, blob_y, blob_radius);
blob_c = repmat(linspace(1500,1600,Ny),Nx,1);
blob_rho = 1000;

% define sound speed and density of medium
medium.sound_speed = c0 * ones(Nx, Ny);          % sound speed [m/s]
medium.density = rho0 * ones(Nx, Ny);            % density [kg/m^3]

assert(size(medium.sound_speed,1) == size(blob_c,1))
assert(size(medium.sound_speed,2) == size(blob_c,2))
medium.sound_speed(blob == 1) = blob_c(blob == 1);

