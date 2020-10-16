% Plane wave imaging example in 2D

clearvars;

% =========================================================================
% SIMULATING THE MEASUREMENTS
% =========================================================================

% create the computational grid
Nx = 256;           % number of grid points in the x (row) direction
Ny = 128;           % number of grid points in the y (column) direction
dx = 40e-3/Nx;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material propertiesÂ
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define a couple of circular scattering regions
scatterer1_radius = 1;          % [grid points]
scatterer1_x = Nx/2;            % [grid points]
scatterer1_y = Ny/2;            % [grid points]
scatterer1_c = 2 * c0;        % sound speed of scatterer [m/s]
scatterer1_rho = 2 * rho0;    % density of scatterer [kg/m^3]
scatterer1 = makeDisc(Nx, Ny, scatterer1_x, scatterer1_y, scatterer1_radius);
% scatterer2_radius = 15;         % [grid points]
% scatterer2_x = Nx/4;            % [grid points]
% scatterer2_y = Ny/2;            % [grid points]
% scatterer2_c = 2 * c0;        % sound speed of scatterer [m/s]
% scatterer2_rho = 2 * rho0;    % density of scatterer [kg/m^3]
% scatterer2 = makeDisc(Nx, Ny, scatterer2_x, scatterer2_y, scatterer2_radius);

% define sound speed of the propagation medium (scatterers random)
medium.sound_speed = c0 * ones(Nx,Ny);          % [m/s]
temp = rand(Nx,Ny);
medium.sound_speed(scatterer1==1) = temp(scatterer1==1) .* scatterer1_c;
% temp = rand(Nx,Ny);
% medium.sound_speed(scatterer2==1) = temp(scatterer2==1) .* scatterer2_c;

% define mass density of the propagation medium (scatterers random)
medium.density = rho0 * ones(Nx,Ny);            % [kg/m^3]
temp = rand(Nx,Ny);
medium.density(scatterer1==1) = temp(scatterer1==1) .* scatterer1_rho;
% temp = rand(Nx,Ny);
% medium.density(scatterer2==1) = temp(scatterer2==1) .* scatterer2_rho;

% create the time array
cfl = 0.2;              % CFL number
t_end = 2*Ny*dy/c0;     % end time of the simulation [s]
kgrid.makeTime(medium.sound_speed,cfl,t_end);
%kgrid.t_array = makeTime(kgrid,medium.sound_speed,cfl,t_end);

% define an apodised planar source 
source.p0 = zeros(Nx, Ny);
sourcewidth = 0.4;
apodisation = getWin(Nx, 'Gaussian','Param', sourcewidth); % see help file for other window options
source.p0(:, pml_size + 1) = apodisation;

% define a sensor array co-aligned with the source
%sensor.mask = source.p0;
sensor_positions = pml_size:(Nx - pml_size);
sensor.mask = zeros(Nx, Ny);
sensor.mask(sensor_positions, pml_size + 1) = 1;

% Calculate the number of sensor elements used
N_sensors = sum(sensor.mask(:));

% run the simulation
inputs = {'PMLSize', pml_size,'PlotLayout', true, 'PlotSim', true};
sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, inputs{:});

% plot the simulated data
% figure
% imagesc(kgrid.t_array*1e6,sensor_positions,sensor_data)
% xlabel('time [\mus]')
% ylabel('sensor number')

% =========================================================================
% IMAGE FORMATION
% =========================================================================

% remove the source from the time series
cut_off = 40;       % time samples to remove from the beginning
sensor_data = [zeros(N_sensors,cut_off) sensor_data(:,cut_off+1:end)];

% window the data
win = getWin(N_sensors, 'Cosine');
sensor_data_apodised = bsxfun(@times, win, sensor_data);
%sensor_data_apodised = sensor_data;

% reconstruct an image using a k-space method
reflection_image = kspaceLineRecon_US(sensor_data_apodised', dx, kgrid.dt, c0);

% % time gain compensation
% tgc_exponent = 80;
% tgc = exp(tgc_exponent * kgrid.t_array' * c0); % exponential
% %tgc = tgc_exponent * kgrid.t_array' * c0; % linear
% reflection_image = bsxfun(@times, tgc, reflection_image);

% envelope detection
reflection_image = transpose(envelopeDetection(reflection_image'));

% log compression
compression_ratio = 3;
reflection_image = logCompression(reflection_image, compression_ratio, true);

% plot the reconstructed image
figure
imagesc((kgrid.t_array * c0/2)*1e3,kgrid.x_vec*1e3,reflection_image')
xlabel('[mm]')
ylabel('[mm]')


