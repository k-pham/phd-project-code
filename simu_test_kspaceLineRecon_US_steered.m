%% set up simulation

% make grid
dx = 10e-6*2;                 % grid point spacing in the x direction [m]
dy = dx;                      % grid point spacing in the y direction [m]
Nx = 1536/2;                  % number of grid points in the x (row) direction
Ny = 1024/2;                  % number of grid points in the y (column) direction

kgrid = kWaveGrid(Nx,dx,Ny,dy);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material properties
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define scatterers
scatterer1_radius = 1;          % [grid points]
scatterer1_x = Nx/2;            % [grid points]
scatterer1_y = Ny/2;            % [grid points]
scatterer1_c = 2 * c0;        % sound speed of scatterer [m/s]
scatterer1_rho = 2 * rho0;    % density of scatterer [kg/m^3]
scatterer1 = makeDisc(Nx, Ny, scatterer1_x, scatterer1_y, scatterer1_radius);

% make medium
medium.sound_speed = c0 * ones(Nx,Ny);
medium.density = rho0 * ones(Nx,Ny);
medium.sound_speed(scatterer1==1) = scatterer1_c;
medium.density(scatterer1==1) = scatterer1_rho;

% create the time array
cfl = 0.2;              % CFL number
t_end = 2*Ny*dy/c0;     % end time of the simulation [s]
kgrid.makeTime(medium.sound_speed,cfl,t_end);

% =========================================================================
% make *angled plane wave* source
% =========================================================================

source_angle    = 5;                        % [deg]
source_centre_x = 0;                        % [m]
source_centre_y = 2.36e-3+kgrid.y_vec(1);   % [m]
source_width_x  = (Nx-2*pml_size)*dx;       % [m]
source_width_y  = source_width_x * tan(deg2rad(source_angle));      % [m]

source_start_point = [source_centre_x-source_width_x/2 , source_centre_y-source_width_y/2];   % [m]
source_end_point   = [source_centre_x+source_width_x/2 , source_centre_y+source_width_y/2];   % [m]

karray = kWaveArray;
karray.addLineElement(source_start_point, source_end_point)
%karray.setArrayPosition([source_centre_x,source_centre_y],source_angle)

source.p_mask = karray.getArrayBinaryMask(kgrid);

amplitude = 10;             % [Pa]

apodisation_width = 0.4;    % proportion of width Nx
apodisation = getWin(Nx, 'Gaussian', 'Param', apodisation_width);

pulse_tpeak      = 20e-9;   % time of pulse peak pressure [s]
pulse_width      = 16e-9;   % FWHM-width of pulse [s]
pulse_variance   = (pulse_width / ( 2*sqrt(2*log(2)) ) )^2;
pulse = gaussian(kgrid.t_array, 1, pulse_tpeak, pulse_variance);

source_signal = amplitude * apodisation * pulse;

source.p = karray.getDistributedSourceSignal(kgrid, source_signal);

% =========================================================================
% =========================================================================

% make sensor
sensor_positions = pml_size:(Nx-pml_size);
sensor.mask = zeros(Nx, Ny);
sensor.mask(sensor_positions, pml_size+1) = 1;


%% run the simulation

inputs = {'PMLSize', pml_size,'PlotLayout', true, 'PlotSim', true};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputs{:});

% plot the simulated data
figure
imagesc(kgrid.t_array*1e6,sensor_positions,sensor_data)
xlabel('time [\mus]')
ylabel('sensor number')





