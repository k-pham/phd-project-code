% Pulse-echo plane-wave US imaging of resolution 27 um wires in a plane 2D

clearvars;

global kgrid t_array dt

file_dir = 'D:\PROJECT\data\simulations\resolution27um\';
beam_positions = {'central', 'offcentre'};

for beam_pos = beam_positions

    beam_pos = beam_pos{1}; %#ok<FXSET>

%% SET UP EXPERIMENT

% create the computational grid
% x_length = 6e-3;            % length of line scan [m]
% y_length = 3e-3;            % depth of field [m]
dx = 2e-6;                  % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
Nx = 2048;    % number of grid points in the x (row) direction
Ny = 768;     % number of grid points in the y (column) direction
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define resolution targets (TUNGSTEN)
num_wires = 7;
num_layers = 1;
wire_radius = 6;                                % [grid points]
wire_xs  = round((1:1:num_wires)*0.5e-3/dx);      % [grid points]
wire_ys  = round((1:1:num_layers)*0.5e-3/dy);     % [grid points]
wire_c   = 5220;                                % sound speed [m/s]
wire_rho = 19300;                               % density [kg/m^3]
wires_layer = zeros(Nx, Ny);

% loop over diff layers of wires
% for idx_layer = 1:num_layers
idx_layer = 1;
    
% make a layer of num_wires wires
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
% source.p0 = zeros(Nx, Ny);
source.p_mask = zeros(Nx, Ny);
source.p_mask(:, pml_size + 1) = 1;
pressure = gaussian(kgrid.t_array,1,20e-9,(5.5e-9)^2);
sourcewidth = 0.4;
apodisation = getWin(Nx+1000, 'Gaussian', 'Param', sourcewidth);
switch beam_pos
    case 'central'
        apodisation = apodisation(501:end-500);
    case 'offcentre'
        apodisation = apodisation(1:end-1000);
end
source.p = apodisation * pressure;
% source.p0(:, pml_size + 1) = apodisation;

% define a sensor array co-aligned with the source
sensor_positions = pml_size:10:(Nx - pml_size);
sensor.mask = zeros(Nx, Ny);
sensor.mask(sensor_positions,pml_size+1) = 1;

% calculate the number of sensor elements used
num_sensors = sum(sensor.mask(:));


%% SIMULATE UP EXPERIMENT

inputs = {'PMLSize', pml_size, 'PlotLayout', true, 'PlotSim', true};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputs{:});

save([file_dir 'layer_' num2str(idx_layer) '_' beam_pos '_sensor_data'], 'sensor_data')

figure
imagesc(kgrid.t_array*1e6,kgrid.x_vec*1e3,sensor_data)
    xlabel('time / \mus')
    ylabel('x position / mm')
    
% end of loop over diff layers of wires
% end


%% IMAGE FORMATION

params.Nx = size(sensor_data,1);
params.Ny = 1;
params.dx = 10*dx;
params.dy = dy;
params.dt = kgrid.dt;

params.trigger_delay        = 0;
params.Nt_zero_pad_source   = 500;
params.Nt_t0_correct        = -262;
params.file_data            = ['111111\resolution27umPlanar_simul_' beam_pos];

reflection_image = reconstruct2dUSimage(sensor_data, params, c0);
% WARNING: kgrid UPDATED

save([file_dir 'layer_' num2str(idx_layer) '_' beam_pos '_reflection_image'], 'reflection_image')

figure
imagesc(kgrid.x_vec*1e3,kgrid.t_array*c0*1e3/2,reflection_image')
axis image


%% peak size analysis

kgrid.makeTime(medium.sound_speed,cfl,t_end);
t_array = kgrid.t_array;
dt = kgrid.dt;

threshold = 70;
peaksInfo = imagePeakFinder(reflection_image, c0, threshold);

save([file_dir 'layer_' num2str(idx_layer) '_' beam_pos '_peaksInfo'], 'peaksInfo')


end









