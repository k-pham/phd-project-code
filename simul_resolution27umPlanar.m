% Pulse-echo plane-wave US imaging of resolution 27 um wires in a plane 2D

clearvars;

% =========================================================================
% SET UP EXPERIMENT
% =========================================================================

% create the computational grid
x_length = 21e-3;           % length of line scan [m]
y_length = 15e-3;           % depth of field [m]
dx = 3e-6;                  % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
Nx = round(x_length/dx);    % number of grid points in the x (row) direction
Ny = round(y_length/dy);    % number of grid points in the y (column) direction
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define resolution targets (TUNGSTEN)
num_wires = 20;
num_layers = 12;
wire_radius = 4;                                % [grid points]
wire_xs  = round((1:1:num_wires)*1e-3/dx);      % [grid points]
wire_ys  = round((1:1:num_layers)*1e-3/dy);     % [grid points]
wire_c   = 5220;                                % sound speed [m/s]
wire_rho = 19300;                               % density [kg/m^3]
wires_layer = zeros(Nx, Ny);

% loop over diff layers of wires
% for idx_layer = 1:num_layers
idx_layer = 1;
    
% make a layer of 20 wires
for idx_wire = 1:num_wires
    wire_x = wire_xs(idx_wire);
    wire_y = wire_ys(idx_layer);
    wires_layer = wires_layer + makeDisc(Nx, Ny, wire_x, wire_y, wire_radius);
end

% define sound speed and density of medium
medium.sound_speed = c0 * ones(Nx, Ny);          % sound speed [m/s]
medium.density = rho0 * ones(Nx, Ny);            % density [kg/m^3]

medium.sound_speed(wires_layer==1) = wire_c;
medium.density(wires_layer==1) = wire_rho;

% create the time array
cfl = 0.2;              % CFL number
t_end = 2*Ny*dy/c0;     % end time of the simulation [s]
kgrid.makeTime(medium.sound_speed,cfl,t_end);

% define a planar apodised (off-central) source
source.p0 = zeros(Nx, Ny);
sourcewidth = 0.4;
apodisation = getWin(Nx, 'Gaussian', 'Param', sourcewidth);
source.p0(:, pml_size + 1) = apodisation;

% define a sensor array co-aligned with the source
sensor_positions = pml_size:(Nx - pml_size);
sensor.mask = zeros(Nx, Ny);
sensor.mask(sensor_positions,pml_size+1) = 1;

% calculate the number of sensor elements used
num_sensors = sum(sensor.mask(:));

% run the simulation
inputs = {'PMLSize', pml_size, 'PlotLayout', true, 'PlotSim', true};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputs{:});

% end of loop over diff layers of wires
% end



