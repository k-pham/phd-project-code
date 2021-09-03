%% ========================================================================
%                               SIMULATION
% =========================================================================

%% set up simulation

% make grid
dx = 10e-6*4;                 % grid point spacing in the x direction [m]
dy = dx;                      % grid point spacing in the y direction [m]
Nx = 1536/4;                  % number of grid points in the x (row) direction
Ny = 1024/4;                  % number of grid points in the y (column) direction

kgrid = kWaveGrid(Nx,dx,Ny,dy);

% define the thickness of the PML [grid points]
pml_size = 20;

% background material properties
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% define scatterers
scatterer1_radius = 10;         % [grid points]
scatterer1_x = Nx/2;            % [grid points]
scatterer1_y = Ny*3/8;          % [grid points]
scatterer1_c = c0;              % sound speed of scatterer [m/s]
scatterer1_rho = 2*rho0;        % density of scatterer [kg/m^3]
scatterer1 = makeDisc(Nx, Ny, scatterer1_x, scatterer1_y, scatterer1_radius);

% make medium
medium.sound_speed = c0 * ones(Nx,Ny);
medium.density = rho0 * ones(Nx,Ny);
medium.sound_speed(scatterer1==1) = scatterer1_c;
medium.density(scatterer1==1) = scatterer1_rho;

% create the time array
cfl = 0.2;                  % CFL number
t_end = 0.6*2*Ny*dy/c0;     % end time of the simulation [s]
kgrid.makeTime(medium.sound_speed,cfl,t_end);

% -------------------------------------------------------------------------
%                   make *angled plane wave* source
% -------------------------------------------------------------------------

source_angle    = 0;                        % [deg]
source_centre_x = 0;                        % [m]
source_centre_y = kgrid.y_vec(1)+pml_size*dy; % 2.36e-3+kgrid.y_vec(1);   % [m]
source_width_x  = (Nx-2*pml_size)*dx;       % [m]
source_width_y  = source_width_x * tan(deg2rad(source_angle));      % [m]

source_start_point = [source_centre_x-source_width_x/2 , source_centre_y-source_width_y/2];   % [m]
source_end_point   = [source_centre_x+source_width_x/2 , source_centre_y+source_width_y/2];   % [m]

karray = kWaveArray;
karray.addLineElement(source_start_point, source_end_point)
%karray.setArrayPosition([source_centre_x,source_centre_y],source_angle)

source.p_mask = karray.getArrayBinaryMask(kgrid);

source_amplitude = 10;             % [Pa]

source_apodisation_width = 0.4;    % proportion of width Nx
source_apodisation = getWin(Nx, 'Gaussian', 'Param', source_apodisation_width);

source_pulse_tpeak      = 20e-9;   % time of pulse peak pressure [s]
source_pulse_width      = 16e-9;   % FWHM-width of pulse [s]
source_pulse_variance   = (source_pulse_width / ( 2*sqrt(2*log(2)) ) )^2;
source_pulse = gaussian(kgrid.t_array, 1, source_pulse_tpeak, source_pulse_variance);

source_signal = source_amplitude * source_apodisation * source_pulse;

source.p = karray.getDistributedSourceSignal(kgrid, source_signal);

% -------------------------------------------------------------------------

% make sensor
sensor_positions = pml_size:(Nx-pml_size);    % sensor points spaced by dx
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

% make sensor kgrid
sensor_kgrid = kWaveGrid(size(sensor_data,1),dx);


%% ========================================================================
%                        PROCESS DATA (T_0 & ANGLE)
% =========================================================================

%% extract TOA from timeseries data

[~,TOA_x_idx] = max(sensor_data,[],2);     % [time index]
TOA_x = TOA_x_idx * kgrid.dt;              % [s]


%% fit plane wave (line) to source

% regression in 1d TOA(x) = p1 + p2*x
Model = [ones(length(sensor_kgrid.x_vec),1) sensor_kgrid.x_vec];
params = Model \ TOA_x;

% evaluate fit over scanning area
TOA_fit = params(1) + params(2) * kgrid.x_vec;

% plot TOA & fit
figure
hold on
plot(sensor_kgrid.x_vec*1e3,TOA_x*1e6,'.')
plot(kgrid.x_vec*1e3,TOA_fit*1e6,'--')
xlabel('x axis / mm')
ylabel('TOA [\mus]')


%% calculate steering angle a
% TOA(x) = p1 + p2*x
% c*p2 = sin(a)

a = asin(c0*params(2));

disp(['steering angle of plane wave is: ' num2str(rad2deg(a)) ' deg'])


%% subtract source from data using saved data

sensor_data_sourceOnly = load('D:\PROJECT\data\simulations\angleComp\test_kspaceLineRecon_US_steered\sensor_data_sourceOnly_40um_0.6t_end.mat','sensor_data_sourceOnly');
sensor_data_sourceOnly = sensor_data_sourceOnly.sensor_data_sourceOnly;

sensor_data = sensor_data - sensor_data_sourceOnly;


%% remove acoustic source from signal by zero padding

% source_pads = 1000;
% sensor_data = cat(2, zeros(kgrid.Ny,source_pads), sensor_data(:,source_pads+1:end));


%% ========================================================================
%                              RECONSTRUCTION
% =========================================================================

%% reconstruction with t0 correction loop

% indeces for central timeseries
nxc = round(sensor_kgrid.Nx/2);

% get TOA(source) in central timeseries
TOA_source = TOA_x_idx(nxc);

% determine true t0 of excitation
t0_excitation_true = source_pulse_tpeak / kgrid.dt;


%% t0 correction loop
% pretend t0 of excitation is unknown and loop through to find out
for t0_excitation = 0%-500:200:500 % [dt]


%% zero padding source offset to sensor plane
% assumes that correction involves zero padding rather than trimming

% set t0(source) in dt
t0_source = t0_excitation - (TOA_source - t0_excitation);

% number of pads
pads = -t0_source;

% check that correction requires zero-padding and not trimming
assert(pads>0)

% apply t0 correction
sensor_data_padded = cat(2, zeros(sensor_kgrid.Nx,pads), sensor_data);


%% reconstruction

% reflection_image = kspaceLineRecon_US_steered(sensor_data_padded,sensor_kgrid.dx,kgrid.dt,c0,a);
reflection_image = kspaceLineRecon_US(sensor_data_padded,sensor_kgrid.dx,kgrid.dt,c0);


%% envelope detection

disp('Envelope detecting ...')
tic
% reflection_image_env = envelopeDetection(reflection_image);
reflection_image_env = transpose(envelopeDetection(reflection_image'));
assert(isequal( size(reflection_image_env), size(reflection_image) ))
disp(['  completed in ' scaleTime(toc)]);


%% plot image

z_vec = kgrid.t_array * c0 / 2;

figure
%imagesc(reflection_image_env')
imagesc(sensor_kgrid.x_vec*1e3,z_vec*1e3,reflection_image_env')
title(['angle ' num2str(source_angle) ', t0 = ' num2str(t0_excitation) '*dt'])
ylabel('depth / mm')
xlabel('x axis / mm')
colorbar
axis image


%%
% pause
end % t0

