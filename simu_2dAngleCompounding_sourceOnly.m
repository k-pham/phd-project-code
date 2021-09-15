%% file directory paths

file_dir_data = 'D:\PROJECT\data\simulations\angleComp\2dAngleCompounding\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\angleComp\2dAngleCompounding\';


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

% % define scatterers
% scatterer1_radius = 5;         % [grid points]
% scatterer1_x = Nx/2;            % [grid points]
% scatterer1_y = Ny*3/8;          % [grid points]
% scatterer1_c = c0;              % sound speed of scatterer [m/s]
% scatterer1_rho = 2*rho0;        % density of scatterer [kg/m^3]
% scatterer1 = makeDisc(Nx, Ny, scatterer1_x, scatterer1_y, scatterer1_radius);

% define scattering medium
scatt_c   = 0;
scatt_rho = 80;
c_rand    = ( rand(Nx,Ny) - 0.5 ) * scatt_c;
rho_rand  = ( rand(Nx,Ny) - 0.5 ) * scatt_rho;

% define non-scattering hole
hole_radius = 15;
hole_x      = Nx/2;
hole_y      = Ny*3/8;
hole_c      = c0;
hole_rho    = rho0;
hole = makeDisc(Nx,Ny,hole_x,hole_y,hole_radius);

% make medium
medium.sound_speed = c0   * ones(Nx,Ny);
medium.density     = rho0 * ones(Nx,Ny);
% medium.sound_speed(scatterer1==1) = scatterer1_c;
% medium.density(scatterer1==1)     = scatterer1_rho;
% medium.sound_speed = medium.sound_speed + c_rand;
% medium.density     = medium.density     + rho_rand;
% medium.sound_speed(hole==1) = hole_c;
% medium.density(hole==1)     = hole_rho;

% create the time array
cfl = 0.2;                  % CFL number
shorten_time = 0.6;         % fraction to shorten length of time array
t_end = shorten_time*2*Ny*dy/c0;     % end time of the simulation [s]
kgrid.makeTime(medium.sound_speed,cfl,t_end);

% -------------------------------------------------------------------------
%                   make *angled plane wave* source
% -------------------------------------------------------------------------
%% anle loop
for source_angle = -30:1:30 %-14:2:14 %[-20:5:20] 
    disp('=================================================')
	disp(['SOURCE ANGLE = ' num2str(source_angle) ' deg'])
% source_angle    = 5;                        % [deg]
source_centre_x = 0;                        % [m]
source_offset_y = 2.36e-3;                  % [m]
source_centre_y = source_offset_y+kgrid.y_vec(1);   % [m] % kgrid.y_vec(1)+pml_size*dy;
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
num_sensors = length(sensor_positions);
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


%% (save) sensor_data for sourceOnly

file_data_sourceOnly = [file_dir_data 'sensor_data_sourceOnly\scattTMM' ...
                        '_' num2str(dx*1e6) 'um' ...
                        '_' num2str(shorten_time) 't_end' ...
                        '_offsetSource' num2str(source_offset_y*1e3) 'e-3' ...
                        '_angle' num2str(source_angle) ...
                        '_sensor_data_sourceOnly.mat'];
sensor_data_sourceOnly = sensor_data;
save(file_data_sourceOnly,'sensor_data_sourceOnly')
end
