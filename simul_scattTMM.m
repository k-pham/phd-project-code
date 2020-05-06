% Pulse-echo plane-wave US imaging - scattering TMM with non-scattering object
% modified to make class-like structures

%% FILE PATHS

% clear all
close all

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';


%% SPECIFY SIMULATION PARAMS -> struct SIMU.PARAMS

% background material properties (WATER)
simu.params.c0   = 1500;    % sound speed [m/s]
simu.params.rho0 = 1000;    % density [kg/m^3]

% choose option to represent scattering medium and object shape
simu.params.scattering_type = 'random';     % options: 'random', 'points'
simu.params.object_shape    = 'hole';       % options: 'hole', 'slab'

% define c/rho for scattering medium and object
simu.params.c_scatt    = 40;        % [m/s]
simu.params.rho_scatt  = 80;        % [kg/m^3]
simu.params.c_object   = 1500;      % [m/s]
simu.params.rho_object = 1000;      % [kg/m^3]

% make medium attenuating (or not)
simu.params.attenuating = true;    % TOGGLE
if simu.params.attenuating
    simu.params.medium_attenuation_coeff = 0.0022;      % [dB MHz^-pow cm^-1]
    simu.params.medium_attenuation_power = 2;           % between 1..3
end

% params for sensor must be set to false here, can change later on
simu.params.sensor_freq_filtered = false;
simu.params.gaussian_freq_filtered = false;
simu.params.sensor_noisy = false;


%% LOAD/MAKE NEW SIMULATION -> struct SIMU

simu_file_path = [file_dir_data file_name(simu.params) '_simu.mat'];

if exist(simu_file_path, 'file')
    disp('Loading existing simulation with specified params..')
    load(simu_file_path, 'simu');
else
    disp('Making new simulation with specified params..')
    simu = make_new_simu(simu);
    save(simu_file_path, 'simu', '-v7.3')
end

fig_simu = plot_simu_medium(simu);
    % saveas(fig_simu,[file_dir_figs file_name(simu) '_medium.fig'])
    % saveas(fig_simu,[file_dir_figs file_name(simu) '_medium.jpg'])


%% OPTION: make *existing* non-attenuating medium attenuating & save simu

if simu.params.attenuating == false
    simu.params.attenuating = true;     % TOGGLE
    simu.params.medium_attenuation_coeff = 0.0022;      % [dB MHz^-pow cm^-1]
    simu.params.medium_attenuation_power = 2;           % between 1..3
    simu = maybe_set_medium_attenuation(simu);
    
    save([file_dir_data file_name(simu.params) '_simu.mat'], 'simu', '-v7.3')
end


%% SPECIFY SIMULATION PARAMS FOR SENSOR -> struct SIMU.PARAMS

% filter with sensor frequency response (or not)
simu.params.sensor_freq_filtered = true;       % TOGGLE

% filter with gaussian frequency filter (or not)
simu.params.gaussian_freq_filtered = false;     % TOGGLE
if simu.params.gaussian_freq_filtered
    simu.params.freq_filter_cf = 1e6;               % [Hz]
    simu.params.freq_filter_bw = 10e6;              % [Hz]
end

% add noise to sensor data
simu.params.sensor_noisy = false;                % TOGGLE


%% LOAD/GENERATE NEW SENSOR DATA -> struct SENSOR

sensor_file_path = [file_dir_data file_name(simu.params) '_sensor.mat'];

if exist(sensor_file_path, 'file')
    disp('Loading existing sensor data with specified params..')
    load(sensor_file_path, 'sensor');
else
    disp('Generating new sensor data with specified params..')
    sensor.data = kspaceFirstOrder2DC(simu.kgrid, simu.medium, simu.source, simu.sensor, simu.inputs{:});
    sensor = set_sensor_params(sensor, simu);
    sensor = make_sensor_kgrid(sensor, simu);
    sensor = make_sensor_tarray(sensor, simu);
    sensor = maybe_sensor_freq_filter(sensor, simu);
    sensor = maybe_gaussian_freq_filter(sensor, simu);
    save(sensor_file_path, 'sensor', '-v7.3')
end

fig_sens = plot_sensor_data(sensor, simu);
    % saveas(fig_sens,  [file_dir_figs file_name(simu) '_data.fig'])
    % saveas(fig_sens,  [file_dir_figs file_name(simu) '_data.jpg'])


%% OPTION: filter *existing* unfiltered sensor data with sensor frequency response & save sensor

if simu.params.sensor_freq_filtered == false
    simu.params.sensor_freq_filtered = true;    % TOGGLE
    sensor = maybe_sensor_freq_filter(sensor, simu);
    
    save([file_dir_data file_name(simu.params) '_sensor.mat'], 'sensor', '-v7.3')
end


%%  OPTION: filter *existing* unfiltered sensor data with gaussian bandpass filter & save sensor

if simu.params.gaussian_freq_filtered == false
    simu.params.gaussian_freq_filtered = true;  % TOGGLE
    simu.params.freq_filter_cf = 1e6;               % [Hz]
    simu.params.freq_filter_bw = 20e6;              % [Hz]
    sensor = maybe_gaussian_freq_filter(sensor, simu);
    
    save([file_dir_data file_name(simu.params) '_sensor.mat'], 'sensor', '-v7.3')
end


%% OPTION: add noise to sensor data before reconstruction

simu.params.sensor_noisy = true;    % TOGGLE


%% RECONSTRUCTION -> struct IMAGE

image.data = reconstruct2dUSimage(sensor.data, sensor.params, simu.params.c0);
    % NOTE: kgrid and t_array UPDATED
    global kgrid t_array

image.kgrid   = kgrid;
image.t_array = t_array;

save([file_dir_data file_name(simu.params) '_image.mat'], 'image' , '-v7.3')
save_image_for_sliceViewer(image, sensor, simu, file_dir_data)

fig_imag = plot_image_data(image, simu);
    % saveas(fig_imag, [file_dir_figs file_name(simu) '_image.fig'])
    % saveas(fig_imag, [file_dir_figs file_name(simu) '_image.jpg'])


%% FUNCTIONS

function pointscatts = get_pointscatt_locations(Nx, Ny, dx, dy, num_points_per_voxel, vox_size)

    num_vox_x = floor(Nx*dx/vox_size);
    num_vox_y = floor(Ny*dy/vox_size);
    vox_x_bounds = round(linspace(1,Nx+1,num_vox_x+1));
    vox_y_bounds = round(linspace(1,Ny+1,num_vox_y+1));

    pointscatts = zeros(Nx,Ny);
    for idx_x = 1 : num_vox_x
        x_left  = vox_x_bounds(idx_x);
        x_right = vox_x_bounds(idx_x+1)-1;
        for idx_y = 1 : num_vox_y
            y_top    = vox_y_bounds(idx_y);
            y_bottom = vox_y_bounds(idx_y+1)-1;

            current_vox_pointscatts = zeros(x_right-x_left+1, y_bottom-y_top+1);
            scatt_locations = randsample(numel(current_vox_pointscatts), num_points_per_voxel);
            current_vox_pointscatts(scatt_locations) = 1;

            pointscatts(x_left:x_right,y_top:y_bottom) = current_vox_pointscatts;

        end
    end

end

function hole = get_hole_location(Nx, Ny)

    hole_radius = 50;               % [grid points]
    hole_x      = round(Nx/2);      % [grid points]
    hole_y      = round(Ny/4);      % [grid points]
    
    hole = makeDisc(Nx,Ny,hole_x,hole_y,hole_radius);
    
end

function slab = get_slab_location(Nx, Ny)

    slab_thickness = 100;           % [grid points]
    slab_y         = 250;           % [grid points]
    slab_yrange    = round(slab_y-slab_thickness/2+0.5 : slab_y+slab_thickness/2-0.5);
    
    slab = zeros(Nx,Ny);
    slab(:,slab_yrange) = 1;
    
end

function filename = file_name(params)
% makes:    file name stem
% requires: simu.params - for scattering medium and object specs

    filename = [params.scattering_type '_SCATT_c'  num2str(params.c_scatt)  '_rho' num2str(params.rho_scatt) '_' ...
                params.object_shape    '_OBJECT_c' num2str(params.c_object) '_rho' num2str(params.rho_object) ];
	
    if params.attenuating
        filename = [filename '_ATTEN_alp' num2str(params.medium_attenuation_coeff) '_pow' num2str(params.medium_attenuation_power) ];
    end
    
    if params.sensor_freq_filtered
        filename = [filename '_FILTER_sfr'];
    end
    
    if params.gaussian_freq_filtered
        filename = [filename '_FILTER_f' num2str(params.freq_filter_cf/1e6) '_bw' num2str(params.freq_filter_bw/1e6) ];
    end

end

function t0_correct = find_t0_correct(sensordata)

    MIP_t   = mean(sensordata(:,1:50),1);
    [~,idx] = max(MIP_t);
    
    t0_correct = -idx +1;

end


%% METHODS FOR SIMU

function simu = make_new_simu(simu)
% makes:    new simu with kgrid, medium, source, sensor, inputs
% requires: simu.params - for medium properties
% NOTE:     order of function calls important due to dependencies

    simu = make_simu_kgrid(simu);
    simu = make_medium(simu);
    simu = maybe_set_medium_attenuation(simu);
    simu = make_time(simu);
    simu = set_pml(simu);
    simu = make_source(simu);
    simu = make_sensor(simu);
    simu = set_inputs(simu);

end

function simu = make_simu_kgrid(simu)
% makes:    simu.kgrid
% requires: nothing

    dx = 10e-6;                 % grid point spacing in the x direction [m]
    dy = dx;                    % grid point spacing in the y direction [m]
    Nx = 1536;                  % number of grid points in the x (row) direction
    Ny = 1024;                  % number of grid points in the y (column) direction
    
    simu.kgrid = kWaveGrid(Nx, dx, Ny, dy);

end

function simu = make_medium(simu)
% makes:    simu.medium.sound_speed/density
% requires: simu.kgrid                  - for size
%           simu.params.c0/rho0         - for background medium
%           simu.params.scattering_type - for type  of scattering medium
%           simu.params.c/rho_scatt     - for c/rho of scattering medium
%           simu.params.object_shape    - for shape of object
%           simu.params.c/rho_object    - for c/rho of object

    % make background medium
    simu.medium.sound_speed = simu.params.c0   * ones(simu.kgrid.Nx, simu.kgrid.Ny);
    simu.medium.density     = simu.params.rho0 * ones(simu.kgrid.Nx, simu.kgrid.Ny);
    
    % make scattering medium specified by scattering_type
    switch simu.params.scattering_type
        case 'random'
            c_rand    = ( rand(simu.kgrid.Nx,simu.kgrid.Ny) - 0.5 ) * simu.params.c_scatt;
            rho_rand  = ( rand(simu.kgrid.Nx,simu.kgrid.Ny) - 0.5 ) * simu.params.rho_scatt;
            
            simu.medium.sound_speed = simu.medium.sound_speed + c_rand;
            simu.medium.density     = simu.medium.density     + rho_rand;
        case 'points'
            num_points_per_voxel = 2;
            vox_size = 100e-6;          % [m]
            pointscatts = get_pointscatt_locations(simu.kgrid.Nx, simu.kgrid.Ny, ...
                                                   simu.kgrid.dx, simu.kgrid.dy, ...
                                                   num_points_per_voxel, vox_size);
            
            simu.medium.sound_speed(pointscatts==1) = simu.params.c_scatt;
            simu.medium.density(pointscatts==1)     = simu.params.rho_scatt;
    end
    
    % make object specified by object_shape
    switch simu.params.object_shape
        case 'hole'
            object = get_hole_location(simu.kgrid.Nx, simu.kgrid.Ny);
        case 'slab'
            object = get_slab_location(simu.kgrid.Nx, simu.kgrid.Ny);
    end
    
    simu.medium.sound_speed(object==1) = simu.params.c_object;
    simu.medium.density(object==1)     = simu.params.rho_object;

end

function simu = maybe_set_medium_attenuation(simu)

    if simu.params.attenuating
        simu.medium.alpha_coeff = simu.params.medium_attenuation_coeff;
        simu.medium.alpha_power = simu.params.medium_attenuation_power;
    end

end

function simu = make_time(simu)
% makes:    simu.kgrid.t_array
% requires: simu.kgrid  - for length of acquisition
%           simu.medium - for max sound speed (to determine dt)

    cfl = 0.2;                  % clf number
    shorten_time = 0.75;        % fraction by which to shorten acquisition length
    t_end = shorten_time * 2 * simu.kgrid.Ny * simu.kgrid.dy / simu.params.c0;      % [s]
    
    simu.kgrid.makeTime(simu.medium.sound_speed, cfl, t_end);

end

function simu = set_pml(simu)
% makes:    simu.params.pml_size
% requires: nothing

    simu.params.pml_size = 20;              % thickness of the PML [grid points]

end

function simu = make_source(simu)
% makes:    simu.source
% requires: simu.kgrid           - for size
%           simu.params.pml_size - for position of source

    amplitude = 10;             % [Pa]
    
    apodisation_width = 0.4;    % proportion of width Nx
    apodisation = getWin(simu.kgrid.Nx, 'Gaussian', 'Param', apodisation_width);
    
    pulse_tpeak      = 20e-9;   % time of pulse peak pressure [s]
    pulse_width      = 16e-9;   % FWHM-width of pulse [s]
    pulse_variance   = (pulse_width / ( 2*sqrt(2*log(2)) ) )^2;
    pulse = gaussian(simu.kgrid.t_array, 1, pulse_tpeak, pulse_variance);
    
    simu.source.p_mask = zeros(simu.kgrid.Nx, simu.kgrid.Ny);
    simu.source.p_mask(:, simu.params.pml_size + 1) = 1;
    
    simu.source.p = amplitude * apodisation * pulse;

end

function simu = make_sensor(simu)
% makes:    simu.sensor
% requires: simu.kgrid           - for size
%           simu.params.pml_size - for position of sensor

    sensor_spacing      = 100e-6;                                   % [m]
    sensor_spacing_grid = round(sensor_spacing / simu.kgrid.dx);    % [grid points]
    sensor_positions_x  = simu.params.pml_size : sensor_spacing_grid : (simu.kgrid.Nx - simu.params.pml_size);
    
    simu.sensor.mask = zeros(simu.kgrid.Nx, simu.kgrid.Ny);
    simu.sensor.mask(sensor_positions_x, simu.params.pml_size+1) = 1;
    
    simu.params.sensor_spacing_grid = sensor_spacing_grid;

end

function simu = set_inputs(simu)
% makes:    simu.inputs
% requires: simu.params.pml_size

    simu.inputs = {'PMLSize', simu.params.pml_size, 'PlotLayout', true, 'PlotSim', true};

end

function fig_simu = plot_simu_medium(simu)
% plots:    simu.medium.sound_speed/density
% requires: simu.medium
%           simu.kgrid  - for axes
%           simu.params - for titles

    fig_simu = figure;
    set(gcf,'Position',[200,200,500,700])
    subplot(2,1,1)
    imagesc(simu.kgrid.x_vec*1e3,simu.kgrid.y_vec*1e3,simu.medium.sound_speed')
        axis image
        title(['sound speed: ' simu.params.scattering_type num2str(simu.params.c_scatt)])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar
    subplot(2,1,2)
    imagesc(simu.kgrid.x_vec*1e3,simu.kgrid.y_vec*1e3,simu.medium.density')
        axis image
        title(['density: ' simu.params.scattering_type num2str(simu.params.rho_scatt)])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar

end


%% METHODS FOR SENSOR

function sensor = set_sensor_params(sensor, simu)
% makes:    sensor.params
% requires: sensor.data - for size
%           simu.params - for sensor spacing
%           simu.kgrid  - for grid size, dt

    sensor.params.Nx = size(sensor.data,1);
    sensor.params.Ny = 1;
    sensor.params.dx = simu.params.sensor_spacing_grid * simu.kgrid.dx;
    sensor.params.dy = simu.kgrid.dy;
    sensor.params.dt = simu.kgrid.dt;
    
    sensor.params.trigger_delay      = 0;
    sensor.params.Nt_zero_pad_source = 50;
    sensor.params.Nt_t0_correct      = -16;
    sensor.params.file_data          = '111111\scattTMM_simul';
    
    assert(sensor.params.Nt_t0_correct == find_t0_correct(sensor.data))

end

function sensor = make_sensor_kgrid(sensor, simu)
% makes:    sensor.kgrid
% requires: sensor.params - for kWaveGrid
%           simu.kgrid    - dt, Nt

    sensor.kgrid    = kWaveGrid(sensor.params.Nx, sensor.params.dx, sensor.params.Ny, sensor.params.dy);
    sensor.kgrid.dt = simu.kgrid.dt;
    sensor.kgrid.Nt = simu.kgrid.Nt;

end

function sensor = make_sensor_tarray(sensor, simu)
% makes:    sensor.t_array
% requires: simu.kgrid.t_array

    sensor.t_array  = simu.kgrid.t_array;

end

function sensor = maybe_sensor_freq_filter(sensor, simu)

    if simu.params.sensor_freq_filtered
        
        % compute the double-sided frequency axis
        Nt = sensor.kgrid.Nt;
        Fs = 1/sensor.kgrid.dt;
        if mod(Nt,2) == 0
            f = (-Nt/2:Nt/2-1) * Fs/Nt;
        else
            f = (-(Nt-1)/2:(Nt-1)/2) * Fs/Nt;
        end
        
        % load sensor frequency response
        data     = load('D:\PROJECT\data\sensorCharac\BK31\191209\sensor_frequency_response.mat');
        f_data   = data.frequency;
        sfr_data = data.sensor_freq_response_mean_norm;
        
        % make sensor frequency response double-sided
        f_data   = [-fliplr(f_data)   f_data(2:end)  ];
        sfr_data = [ fliplr(sfr_data) sfr_data(2:end)];
        
        % resample sfr at f
        sfr = interp1(f_data, sfr_data, f);
        
        % zeropad above 80 MHz to remove noisy sfr & NaNs from resampling
        sfr(abs(f)>80e6) = 0;
        
        % apply filter with zero-phase
        sensor.data = real(ifft(ifftshift( ...
                        bsxfun(@times, sfr, ...
                            fftshift(fft(sensor.data, [], 2), 2) ...
                        ) ...
                      , 2), [], 2));
    end

end

function sensor = maybe_gaussian_freq_filter(sensor, simu)

    if simu.params.gaussian_freq_filtered
        cf = simu.params.freq_filter_cf;
        bw = simu.params.freq_filter_bw;
        
        disp(['FILTER ' num2str(cf/1e6) ' bw ' num2str(bw/1e6)])
        
        win = getWin(sensor.kgrid.Nt,'Tukey');
        sensor.data = sensor.data .* win';
        
        bandwidth_pc = bw / cf * 100;
        sensor.data = gaussianFilter(sensor.data, 1/sensor.kgrid.dt, cf, bandwidth_pc, false);
    end

end

% TO DO: make function to introduce instrument noise

function fig_sens = plot_sensor_data(sensor, simu)
% plots:    sensor.data
% requires: sensor.data
%           sensor.kgrid   - for axes
%           sensor.t_array - for axes
%           simu.params    - for titles

    % figure
    % imagesc(sensor.data(:,1:50)')
    %     xlabel('x position / dx')
    %     ylabel('time / dt')
    
    fig_sens = figure;
    imagesc(sensor.kgrid.x_vec*1e3,sensor.t_array(50:end)*1e6,sensor.data(:,50:end)')
        title([simu.params.scattering_type ' c ' num2str(simu.params.c_scatt) ' rho ' num2str(simu.params.rho_scatt)])
        xlabel('x position / mm')
        ylabel('time / \mus')
        colorbar

end


%% METHODS FOR IMAGE

function fig_imag = plot_image_data(image, simu)
% plots:    image.data
% requires: image.data
%           image.kgrid   - for axes
%           image.t_array - for axes
%           simu.params   - for c0, titles

    fig_imag = figure;
    imagesc(image.kgrid.x_vec*1e3, image.t_array*simu.params.c0*1e3, image.data')
                     % omit factor 1/2 in depth because of doubled depth bug
        axis image
        title([simu.params.scattering_type ' c ' num2str(simu.params.c_scatt) ' rho ' num2str(simu.params.rho_scatt)])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar

end

function save_image_for_sliceViewer(image, sensor, simu, file_dir_data)
% saves:    image.data reshaped to 3D
% requires: image.data
%           image.kgrid   - for volume spacing
%           sensor.params - for dt
%           simu.params   - for c0, file name

    [Nx_image, Ny_image] = size(image.data);
    volume_data = reshape(image.data, Nx_image, 1, Ny_image);
    volume_spacing = [image.kgrid.dx, 1, sensor.params.dt*simu.params.c0];
                     % omit factor 1/2 in dz because of doubled depth bug
    
    save([file_dir_data file_name(simu) '_image_4sliceViewer.mat'], 'volume_data', 'volume_spacing', '-v7.3')

end

