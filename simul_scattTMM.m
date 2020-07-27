% Pulse-echo plane-wave US imaging - scattering TMM with non-scattering object
% modified to make class-like structures

%% (0) FILE PATHS

% clear all
close all

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';


%% (1-SPECIFY): simulation params -> struct SIMU.PARAMS

% background material properties (WATER)
simu.params.c0   = 1500;    % sound speed [m/s]
simu.params.rho0 = 1000;    % density [kg/m^3]

% define scattering medium
simu.params.scatt_type = 'random';          % options: 'random', 'points', 'non-scattering'
simu.params.scatt_c    = 40;                         % [m/s]
simu.params.scatt_rho  = 80;                         % [kg/m^3]

% define object
simu.params.object_shape = 'hole';             % options: 'hole', 'slab', 'wire', 'no object'
simu.params.object_c     = 1500;                    % [m/s]
simu.params.object_rho   = 1000;                    % [kg/m^3]
simu.params.object_x     = 1536/2;                  % [grid points]
simu.params.object_y     = 1024/4;                  % [grid points]

% make medium attenuating (or not)
simu.params.attenuating = false;    % TOGGLE
if simu.params.attenuating
    simu.params.medium_attenuation_coeff = 0.15;      % [dB MHz^-pow cm^-1]
    simu.params.medium_attenuation_power = 1.21;           % between 1..3
end

% shorten time of acquisition
simu.params.shorten_time = 1;                       % [fraction]

% sensor spacing
simu.params.sensor_spacing = 100e-6;                % [m]

% params for sensor must be set to false here, can change later on
simu.params.sensor_freq_filtered = false;
simu.params.gaussian_freq_filtered = false;
simu.params.custom_freq_filtered = false;
simu.params.sensor_noisy = false;

% params for recon must be set to false here, can change later on
simu.params.freq_compound = false;


%% (1) LOAD/MAKE NEW SIMULATION -> struct SIMU

simu_file_path = [file_dir_data file_name(simu.params) '_simu.mat'];

if exist(simu_file_path, 'file')
    disp(['Loading existing simulation: ' file_name(simu.params)])
    load(simu_file_path, 'simu');
else
    disp(['Making new simulation: ' file_name(simu.params)])
    simu = make_new_simu(simu);
    save(simu_file_path, 'simu', '-v7.3')

    fig_simu = plot_simu_medium(simu);
        saveas(fig_simu,[file_dir_figs file_name(simu.params) '_medium.fig'])
        saveas(fig_simu,[file_dir_figs file_name(simu.params) '_medium.jpg'])
end


%% (1-OPTION): make *existing* non-attenuating medium attenuating & save simu
% 
% if simu.params.attenuating == false
%     simu.params.attenuating = true;     % TOGGLE
%     simu.params.medium_attenuation_coeff = 0.9;      % [dB MHz^-pow cm^-1]
%     simu.params.medium_attenuation_power = 1.1;           % between 1..3
%     simu = maybe_set_medium_attenuation(simu);
%     
%     save([file_dir_data file_name(simu.params) '_simu.mat'], 'simu', '-v7.3')
% end


%% (2-SPECIFY): simulation params for sensor -> struct SIMU.PARAMS

% filter with sensor frequency response (or not)
simu.params.sensor_freq_filtered = false;       % TOGGLE

% % filter with gaussian frequency filter (or not)
% simu.params.gaussian_freq_filtered = false;     % TOGGLE
% if simu.params.gaussian_freq_filtered
%     simu.params.freq_filter_cf = 1e6;               % [Hz]
%     simu.params.freq_filter_bw = 20e6;              % [Hz]
% end

% % filter with custom frequency filter (or not)
% simu.params.custom_freq_filtered = false;       % TOGGLE
% if simu.params.custom_freq_filtered
%     % REQUIRES: manual definition of custom filter parameters
%     simu.params.custom_freq_filter_data = {f_custom, filter_custom};	% [Hz]
% end

% add noise to sensor data (or not)
simu.params.sensor_noisy = false;                % TOGGLE
if simu.params.sensor_noisy
    simu.params.sensor_snr = 20;                    % [dB w.r.t. rms]
end


%% (2) LOAD/GENERATE NEW SENSOR DATA -> struct SENSOR

sensor_file_path = [file_dir_data file_name(simu.params) '_sensor.mat'];

if exist(sensor_file_path, 'file')
    disp(['Loading existing sensor data: ' file_name(simu.params)])
    load(sensor_file_path, 'sensor');
else
    disp(['Generating new sensor data: ' file_name(simu.params)])
    sensor = generate_new_sensordata(simu);
    save(sensor_file_path, 'sensor', '-v7.3')
    
    fig_sens = plot_sensor_data(sensor, simu);
        saveas(fig_sens, [file_dir_figs file_name(simu.params) '_data.fig'])
        saveas(fig_sens, [file_dir_figs file_name(simu.params) '_data.jpg'])
end


%% (2) subtract source contribution/background from sensor data -> struct SENSOR

disp('Subtract source contribution from sensor.data')

background_path = ['D:\PROJECT\data\simulations\scattTMM\non-scattering _no object_\' ...
                     num2str(simu.params.sensor_spacing*1e6) ' um\' ...
                     'non-scattering_SCATT_c0_rho0_no object_OBJECT_c0_rho0_x0_y0'];

background_sensor = load([background_path '_sensor.mat']);
bckgr_sensor      = background_sensor.sensor;

background_simu   = load([background_path '_simu.mat']);
bckgr_simu        = background_simu.simu;

% check sensor spacing the same
assert(simu.params.sensor_spacing == bckgr_simu.params.sensor_spacing)

% resample background data if necessary
if simu.kgrid.dt ~= bckgr_simu.kgrid.dt
    bckgr_sensor.data = permute(interp1(bckgr_sensor.t_array, bckgr_sensor.data', sensor.t_array, 'pchip'), [2 1]);
end

% subtract background data
sensor_data = sensor.data;      % save background-unsubtracted data for 2dfft
sensor.data = sensor.data - bckgr_sensor.data;

file_data_bckgrsubtr = [file_dir_figs file_name(simu.params) '_data_bckgrsubtracted'];
if ~exist([file_data_bckgrsubtr '.fig'],'file')
fig_sens2 = plot_sensor_data(sensor, simu);
    saveas(fig_sens2, [file_data_bckgrsubtr '.fig'])
    saveas(fig_sens2, [file_data_bckgrsubtr '.jpg'])
end


%% (2-OPTION): filter *existing* unfiltered sensor data with sensor frequency response & save sensor
% 
% if simu.params.sensor_freq_filtered == false
%     simu.params.sensor_freq_filtered = true;    % TOGGLE
%     sensor = maybe_sensor_freq_filter(sensor, simu);
%     
%     save([file_dir_data file_name(simu.params) '_sensor.mat'], 'sensor', '-v7.3')
% end


%% MIGRATED TO 3-SPECIFY ... (2-OPTION): filter *existing* unfiltered sensor data with gaussian bandpass filter & save sensor
% % 
% % if simu.params.gaussian_freq_filtered == false
% %     simu.params.gaussian_freq_filtered = true;  % TOGGLE
% %     simu.params.freq_filter_cf = 1e6;               % [Hz]
% %     simu.params.freq_filter_bw = 20e6;              % [Hz]
% %     sensor = maybe_gaussian_freq_filter(sensor, simu);
% %     
% %     save([file_dir_data file_name(simu.params) '_sensor.mat'], 'sensor', '-v7.3')
% % end


%% MIGRATED TO 3-SPECIFY ... (2-OPTION): filter *existing* unfiltered sensor data with custom filter - here: compound bandpass filter
% % 
% % if simu.params.custom_freq_filtered == false
% %     simu.params.custom_freq_filtered = true;    % TOGGLE
% %     
% %     % manually define custom filter parameters for compound bandpass filter:
% %     centre_freqs = (1:1:10)*1e6;        % [Hz]
% %     bandwidth    = 2e6;                 % [Hz]
% %     [f_compound, filter_compound] = get_compound_filter(centre_freqs, bandwidth);
% %     
% %     simu.params.custom_freq_filter_cfs  = centre_freqs;
% %     simu.params.custom_freq_filter_bw   = bandwidth;
% %     simu.params.custom_freq_filter_data = {f_compound, filter_compound};	% [Hz]
% %     % sensor = maybe_custom_freq_filter(sensor, simu);
% %     
% %     % save([file_dir_data file_name(simu.params) '_sensor.mat'], 'sensor', '-v7.3')
% % end


%% (2-OPTION): add noise to *existing* non-noisy sensor data before reconstruction & save sensor
% 
% if simu.params.sensor_noisy == false
%     simu.params.sensor_noisy = true;    % TOGGLE
%     simu.params.sensor_snr = 20;                    % [dB w.r.t. rms]
%     sensor = maybe_make_sensor_noisy(sensor, simu);
%     
%     save([file_dir_data file_name(simu.params) '_sensor.mat'], 'sensor', '-v7.3')
% end


%% (3-SPECIFY): simulation params for reconstruction -> struct SIMU.PARAMS

% recon with gaussian frequency filter
simu.params.gaussian_freq_filtered = false;
if simu.params.gaussian_freq_filtered
    simu.params.freq_filter_cf = 1e6;               % [Hz]
    simu.params.freq_filter_bw = 20e6;              % [Hz]
end

% recon with custom (compound) frequency filter
simu.params.custom_freq_filtered = false;
if simu.params.custom_freq_filtered
    % manually define custom filter parameters for compound bandpass filter:
    centre_freqs = (1:1:10)*1e6;        % [Hz]
    bandwidth    = 2e6;                 % [Hz]
    [f_compound, filter_compound] = get_compound_filter(centre_freqs, bandwidth);
    
    simu.params.custom_freq_filter_cfs  = centre_freqs;
    simu.params.custom_freq_filter_bw   = bandwidth;
    simu.params.custom_freq_filter_data = {f_compound, filter_compound};    % [Hz]
end

% recon with frequency compounding
simu.params.freq_compound = false;
if simu.params.freq_compound
    simu.params.freq_compound_method = 'incoherent';      % options: 'coherent', 'incoherent'
    simu.params.freq_compound_cf     = (1:1:10)*1e6;    % [Hz]
    simu.params.freq_compound_bw     = 2e6;             % [Hz]
end


%% (3) LOAD/RECONSTRUCT NEW IMAGE -> struct IMAGE

image_file_path = [file_dir_data file_name(simu.params) '_image.mat'];

if exist(image_file_path, 'file')
    disp(['Loading existing image data: ' file_name(simu.params)])
    load(image_file_path, 'image');
else
    disp(['Reconstructing new image data: ' file_name(simu.params)])
    image = recon_new_image(sensor, simu);
    save(image_file_path, 'image', '-v7.3')
    % save_image_for_sliceViewer(image, sensor, simu, [file_dir_data file_name(simu.params)])
    
    fig_imag = plot_image_data(image, simu);
        saveas(fig_imag, [file_dir_figs file_name(simu.params) '_image.fig'])
        saveas(fig_imag, [file_dir_figs file_name(simu.params) '_image.jpg'])
end

image.c0 = simu.params.c0;


%% (4) ANALYSE IMAGE

switch simu.params.object_shape
    
    % scattering distributions in hole & stmm, plot if wanted, scattering SNR / CNR
    case 'hole'
        
        plot_toggle = true;
        
        if plot_toggle == true
            fig_distr = figure;
                hold on
                xlabel('pixel intensity')
                ylabel('count')
        end
        
        [scatter_hole_mean, scatter_hole_std] = get_scattering_distr_in_hole(image, simu, plot_toggle);
        [scatter_stmm_mean, scatter_stmm_std] = get_scattering_distr_in_stmm(image, simu, plot_toggle);
        
        scatSNR = scatter_stmm_mean / scatter_hole_mean;
        scatCNR = (scatter_stmm_mean - scatter_hole_mean) / (scatter_stmm_std + scatter_hole_std);
        
        if plot_toggle == true
            title(['scattering distributions, SNR = ' num2str(scatSNR) ', CNR = ' num2str(scatCNR)])
            legend(gca,'show')
            saveas(fig_distr, [file_dir_figs file_name(simu.params) '_image_distr.fig'])
            saveas(fig_distr, [file_dir_figs file_name(simu.params) '_image_distr.jpg'])
        end
        
        disp('  scatSNR   scatCNR')
        disp([scatSNR,scatCNR])
        
    % specular signal & resolution
    case 'wire'
        
        signal_wire = get_specular_signal_of_wire(image, simu);
        [resoLat, resoAxi] = get_resolution_of_wire(image, simu);
        
        disp('  specSig   resoLat   resoAxi')
        disp([signal_wire,resoLat*1e6,resoAxi*1e6])
        
end


%% (4) ANALYSE SENSOR DATA: frequency content in 2d-fft

disp('Fourier transforming data in 2D..')

[size_x, size_t] = size(sensor_data);       % use background-unsubtracted data and subtract 2dfft background data later

size_t_fft = round((size_t+1)/2);
size_x_fft = round((size_x+1)/2);

sensor_data_fftT  = zeros(size_x    , size_t_fft);
sensor_data_fftTX = zeros(size_x_fft, size_t_fft);

for x = 1 : sensor.params.Nx
    [freqT, sensor_data_fftT(x,:) ] = spect(sensor_data(x,:), 1/sensor.params.dt);
end

for t = 1 : size_t_fft
    [freqX, sensor_data_fftTX(:,t)] = spect(sensor_data_fftT(:,t), 1/sensor.params.dx);
end

% cut out super high freq
omegarange = 1:round(length(freqT)/2);
freqT = freqT(omegarange);
sensor_data_fftTX = sensor_data_fftTX(:,omegarange);

% file location for background 2dfft data
file_background_2dfft = ['D:\PROJECT\data\simulations\scattTMM\non-scattering _no object_\' ...
                            num2str(simu.params.sensor_spacing*1e6) ' um\sensor_data_2dfft.mat'];

% if empty target: save background 2dfft data
if strcmp(simu.params.scatt_type,'non-scattering') && strcmp(simu.params.object_shape,'no object')
    % only save if not already exists
    if ~exist(file_background_2dfft,'file')
        disp('Saving background 2dfft data..')
        sensor_data_fftTX_background = sensor_data_fftTX;
        freqT_background = freqT;
        freqX_background = freqX;
        save(file_background_2dfft, 'sensor_data_fftTX_background', 'freqT_background', 'freqX_background')
    else
        disp('Background 2dfft data already exists.')
    end
    
% else: subtract it from current data
else
    disp('Subtracting background 2dfft data..')
    
    % loading bacground 2dfft data
    load(file_background_2dfft, 'sensor_data_fftTX_background', 'freqT_background', 'freqX_background');
    
    % resample background 2dfft data onto same grid as current 2dfft data
    sensor_data_fftTX_background_resample = permute(interp1(freqT_background, sensor_data_fftTX_background', freqT), [2 1]);
    
    % subtract background 2dfft data
    sensor_data_fftTX = sensor_data_fftTX - sensor_data_fftTX_background_resample;
    
    % plot 2dfft data
    x_min = 4;
    fig_2dfft = figure('Position',[300,300,750,450]);
    imagesc(freqT/1e6, freqX(x_min:end)/1e3, sensor_data_fftTX(x_min:end,:))
        title('2D FFT of sensor.data')
        xlabel('Temporal frequency \omega [MHz]')
        ylabel('Spatial frequency k_x [mm^{-1}]')
        xlim([0,70])
    	ylim([0,20])
        colorbar
        set(gca,'FontSize',13)
        % caxis([-0.5e-5, 1e-5])
    saveas(fig_2dfft, [file_dir_figs file_name(simu.params) '_sensor_2dfft.fig'])
	saveas(fig_2dfft, [file_dir_figs file_name(simu.params) '_sensor_2dfft.jpg'])
    
end


%% transform 2dfft data to angle of incidence \theta dependence

disp('Transforming 2dfft data to angle-of-incidence dependence..')

c = 1500;

% make sin(theta) matrix and make non-physical values NaNs
sintheta = freqX' * c ./ freqT;
sintheta(sintheta>1) = NaN;

% convert to theta in degrees
theta = asin(sintheta);
theta = theta / pi * 180;

% % plot theta matrix in k-space axes
% figure('Position',[300,300,750,450]);
% imagesc(freqT/1e6, freqX/1e3, theta)
%     title('theta')
%     xlabel('Temporal frequency \omega [MHz]')
%     ylabel('Spatial frequency k_x [mm^{-1}]')
%     xlim([0,70])
% 	ylim([0,5])
%     colorbar
%     set(gca,'FontSize',13)

% mask of theta is not non-physical/NaN
mask = not(isnan(theta));

% extract theta, freqT and sensor_data_fftTX at mask
theta_ar             = theta(mask);
sensor_data_fftTX_ar = sensor_data_fftTX(mask);
freqT_mat            = repmat(freqT, [length(freqX), 1]);
freqT_ar             = freqT_mat(mask);

% % plot freqT matrix in k-space axes
% figure('Position',[300,300,750,450]);
% imagesc(freqT/1e6, freqX/1e3, freqT_mat/1e6)
%     title('freqT \omega [MHz]')
%     xlabel('Temporal frequency \omega [MHz]')
%     ylabel('Spatial frequency k_x [mm^{-1}]')
%     xlim([0,70])
%     ylim([0,5])
%     caxis([0,70])
%     colorbar
%     set(gca,'FontSize',13)

% make function for interpolating sensor_data_fftTX
F = scatteredInterpolant(theta_ar, freqT_ar, sensor_data_fftTX_ar);

% resample sensor_data_fftTX at regular theta and freqT arrays
theta_new = 0:0.1:90;
freqT_new = (0:0.1:75)*1e6;
sensor_data_fftTA = F({theta_new, freqT_new});

% % plot sensor_data_fft(theta,freqT) as scatter plot
% figure
% scatter(freqT_ar/1e6, theta_ar, 1, sensor_data_fftTX_ar)
%     title('sensor.data')
%     xlabel('Temporal frequency \omega [MHz]')
%     ylabel('Angle of incidence [deg]')
%     xlim([0,75])
%     ylim([0,90])
%     colorbar
%     set(gca,'FontSize',13)

% plot sensor_data_fft(theta,freqT) as interpolated scatter data
f_min = 5;
fig_2dfft_aoi = figure('Position',[300,300,750,450]);
imagesc(freqT_new(f_min:end)/1e6, theta_new, sensor_data_fftTA(:,f_min:end))
    title('sensor.data')
    xlabel('Temporal frequency \omega [MHz]')
    ylabel('Angle of incidence [deg]')
    xlim([0,75])
    ylim([0,90])
    colorbar
    set(gca,'FontSize',13)
saveas(fig_2dfft_aoi, [file_dir_figs file_name(simu.params) '_sensor_2dfft_aoi.fig'])
saveas(fig_2dfft_aoi, [file_dir_figs file_name(simu.params) '_sensor_2dfft_aoi.jpg'])


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

function hole = get_hole_location(Nx, Ny, hole_x, hole_y)

    hole_radius = 50;               % [grid points]
    
    hole = makeDisc(Nx,Ny,hole_x,hole_y,hole_radius);
    
end

function slab = get_slab_location(Nx, Ny)

    slab_thickness = 100;           % [grid points]
    slab_y         = 250;           % [grid points]
    slab_yrange    = round(slab_y-slab_thickness/2+0.5 : slab_y+slab_thickness/2-0.5);
    
    slab = zeros(Nx,Ny);
    slab(:,slab_yrange) = 1;
    
end

function wire = get_wire_location(Nx, Ny, wire_x, wire_y)

    wire_radius = 0;                % [grid points]
    
    wire = makeDisc(Nx,Ny,wire_x,wire_y,wire_radius);

end

function filename = file_name(params)
% makes:    file name stem
% requires: simu.params - for scattering medium and object specs

    filename = [params.scatt_type    '_SCATT_c'  num2str(params.scatt_c)  '_rho' num2str(params.scatt_rho) '_' ...
                params.object_shape  '_OBJECT_c' num2str(params.object_c) '_rho' num2str(params.object_rho) '_' ...
                                             'x' num2str(params.object_x) '_y'   num2str(params.object_y) ];
	
    if params.attenuating
        filename = [filename '_ATTEN_alp' num2str(params.medium_attenuation_coeff) '_pow' num2str(params.medium_attenuation_power) ];
    end
    
    if params.sensor_freq_filtered
        filename = [filename '_FILTER_sfr'];
    end
    
    if params.gaussian_freq_filtered
        filename = [filename '_FILTER_f' num2str(params.freq_filter_cf/1e6) '_bw' num2str(params.freq_filter_bw/1e6) ];
    end
    
    if params.custom_freq_filtered
        filename = [filename '_CFILTER_f' num2str(params.custom_freq_filter_cfs(1)/1e6) '-' num2str(params.custom_freq_filter_cfs(end)/1e6) '_bw' num2str(params.custom_freq_filter_bw/1e6) ];
    end
    
    if params.sensor_noisy
        filename = [filename '_NOISE_snr' num2str(params.sensor_snr) ];
    end
    
    if params.freq_compound
        filename = [filename '_COMPOUND_' params.freq_compound_method '_f' num2str(params.freq_compound_cf(1)/1e6) '-' num2str(params.freq_compound_cf(end)/1e6) '_bw' num2str(params.freq_compound_bw/1e6) ];
    end

end

function t0_correct = find_t0_correct(sensordata)

    MIP_t   = mean(sensordata(:,1:50),1);
    [~,idx] = max(MIP_t);
    
    t0_correct = -idx +1;

end

function [f, compound_filter] = get_compound_filter(centre_freqs, bandwidth)

    f         = (0:0.05:400)*1e6;     % [Hz]
    magnitude = 1;
    variance  = (bandwidth / (2 * sqrt(2 * log(2)))).^2;

    figure, hold on
    for mean = centre_freqs
        gauss_filter = max(gaussian(f, magnitude, mean, variance), gaussian(f, magnitude, -mean, variance));

        plot(f/1e6, gauss_filter, ':')
            xlabel('Frequency / MHz')

        if ~exist('compound_filter','var')
            compound_filter = gauss_filter;
        else
            compound_filter = compound_filter + gauss_filter;
        end
    end

    compound_filter = compound_filter/max(compound_filter);
    
    plot(f/1e6, compound_filter, 'k-')
    title('compound filter')
    xlim([0,30])
    
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
%           simu.params.scatt_type      - for type  of scattering medium
%           simu.params.scatt_c/rho     - for c/rho of scattering medium
%           simu.params.object_shape    - for shape of object
%           simu.params.object_c/rho    - for c/rho of object

    % make background medium
    simu.medium.sound_speed = simu.params.c0   * ones(simu.kgrid.Nx, simu.kgrid.Ny);
    simu.medium.density     = simu.params.rho0 * ones(simu.kgrid.Nx, simu.kgrid.Ny);
    
    % make scattering medium specified by scatt_type
    switch simu.params.scatt_type
        case 'random'
            c_rand    = ( rand(simu.kgrid.Nx,simu.kgrid.Ny) - 0.5 ) * simu.params.scatt_c;
            rho_rand  = ( rand(simu.kgrid.Nx,simu.kgrid.Ny) - 0.5 ) * simu.params.scatt_rho;
            
            simu.medium.sound_speed = simu.medium.sound_speed + c_rand;
            simu.medium.density     = simu.medium.density     + rho_rand;
        case 'points'
            num_points_per_voxel = 2;
            vox_size = 100e-6;          % [m]
            pointscatts = get_pointscatt_locations(simu.kgrid.Nx, simu.kgrid.Ny, ...
                                                   simu.kgrid.dx, simu.kgrid.dy, ...
                                                   num_points_per_voxel, vox_size);
            
            simu.medium.sound_speed(pointscatts==1) = simu.params.scatt_c;
            simu.medium.density(pointscatts==1)     = simu.params.scatt_rho;
        case 'non-scattering'
            % do nothing, keep uniform background medium properties
    end
    
    % make object specified by object_shape
    switch simu.params.object_shape
        case 'hole'
            object = get_hole_location(simu.kgrid.Nx, simu.kgrid.Ny, simu.params.object_x, simu.params.object_y);
        case 'slab'
            object = get_slab_location(simu.kgrid.Nx, simu.kgrid.Ny);
        case 'wire'
            object = get_wire_location(simu.kgrid.Nx, simu.kgrid.Ny, simu.params.object_x, simu.params.object_y);
        case 'no object'
            object = zeros(simu.kgrid.Nx, simu.kgrid.Ny);
    end
    
    simu.medium.sound_speed(object==1) = simu.params.object_c;
    simu.medium.density(object==1)     = simu.params.object_rho;

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
    t_end = simu.params.shorten_time * 2 * simu.kgrid.Ny * simu.kgrid.dy / simu.params.c0;      % [s]
    
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

    sensor_spacing_grid = round(simu.params.sensor_spacing / simu.kgrid.dx);    % [grid points]
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
        title(['sound speed: ' simu.params.scatt_type   num2str(simu.params.scatt_c) ', ' ...
                               simu.params.object_shape num2str(simu.params.object_c) ])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar
    subplot(2,1,2)
    imagesc(simu.kgrid.x_vec*1e3,simu.kgrid.y_vec*1e3,simu.medium.density')
        axis image
        title(['density: ' simu.params.scatt_type   num2str(simu.params.scatt_rho) ', ' ...
                           simu.params.object_shape num2str(simu.params.object_rho) ])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar
    
end


%% METHODS FOR SENSOR

function sensor = generate_new_sensordata(simu)
% makes:    new sensor with data, params, kgrid, t_array
% requires: simu
% NOTE:     order of function calls important due to dependencies

    sensor.data = kspaceFirstOrder2DC(simu.kgrid, simu.medium, simu.source, simu.sensor, simu.inputs{:});
    sensor = set_sensor_params(sensor, simu);
    sensor = make_sensor_kgrid(sensor, simu);
    sensor = make_sensor_tarray(sensor, simu);
    sensor = maybe_sensor_freq_filter(sensor, simu);
    sensor = maybe_gaussian_freq_filter(sensor, simu);
    sensor = maybe_custom_freq_filter(sensor, simu);
    sensor = maybe_make_sensor_noisy(sensor, simu);

end

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
    sensor.params.Nt_zero_pad_source = 100;
    sensor.params.Nt_t0_correct      = find_t0_correct(sensor.data);  % hole -16
    sensor.params.file_data          = '111111\scattTMM_simul';
    
    % assert(sensor.params.Nt_t0_correct == find_t0_correct(sensor.data))

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

function sensor = maybe_custom_freq_filter(sensor, simu)

    if simu.params.custom_freq_filtered
        
        % compute the double-sided frequency axis
        Nt = sensor.kgrid.Nt;
        Fs = 1/sensor.kgrid.dt;
        if mod(Nt,2) == 0
            f = (-Nt/2:Nt/2-1) * Fs/Nt;
        else
            f = (-(Nt-1)/2:(Nt-1)/2) * Fs/Nt;
        end
        
        % load custom filter from simu.params
        [f_custom, filter_custom] = simu.params.custom_freq_filter_data{:};
        
        % make custom filter double-sided
        f_custom   = [-fliplr(f_custom)   f_custom(2:end)  ];
        filter_custom = [ fliplr(filter_custom) filter_custom(2:end)];
        
        % resample custom filter at f
        filter = interp1(f_custom, filter_custom, f);
        
        % apply filter with zero-phase
        sensor.data = real(ifft(ifftshift( ...
                        bsxfun(@times, filter, ...
                            fftshift(fft(sensor.data, [], 2), 2) ...
                        ) ...
                      , 2), [], 2));
    end

end

function sensor = maybe_make_sensor_noisy(sensor, simu)

    if simu.params.sensor_noisy
        
        % save source for later addition
        sensor_data_source = sensor.data(:,1:50);
        
        % zero-pad source (so that addNoise can determine signal rms)
        sensor.data(:,1:50) = 0;
        
        % add noise
        sensor.data = addNoise(sensor.data, simu.params.sensor_snr, 'rms');
        
        % add source back in
        sensor.data(:,1:50) = sensor.data(:,1:50) + sensor_data_source;
        
    end
    
end

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
        title([simu.params.scatt_type   ' c ' num2str(simu.params.scatt_c)  ' rho ' num2str(simu.params.scatt_rho) ', ' ...
               simu.params.object_shape ' c ' num2str(simu.params.object_c) ' rho ' num2str(simu.params.object_rho) ])
        xlabel('x position / mm')
        ylabel('time / \mus')
        colorbar
    
end


%% METHODS FOR IMAGE

function image = recon_new_image(sensor, simu)

    % recon with frequency compounding
    if simu.params.freq_compound
        
        method = simu.params.freq_compound_method;
        cfs    = simu.params.freq_compound_cf;
        bw     = simu.params.freq_compound_bw;
        
        % set envelope detection according to compounding method
        switch method
            case 'coherent'
                toEnvelopeDetect = false;
            case 'incoherent'
                toEnvelopeDetect = true;
        end
        
        % loop over centre frequencies and compound
        num_images_compounded = 0;
        for cf = cfs
            disp(['CENTRE FREQ ' num2str(cf/1e6) ', BANDWIDTH ' num2str(bw/1e6)])
            single_image = reconstruct2dUSimage(sensor.data, sensor.params, simu.params.c0, ...
                                          'Upsample', true, ...
                                          'FreqBandFilter', {cf, bw}, ...
                                          'EnvelopeDetect', toEnvelopeDetect, ...
                                          'SaveImageToFile', false );
            % add to compound or create new if first image
            if ~exist('compound_image','var')
                compound_image = single_image;
            else
                compound_image = compound_image + single_image;
            end
            num_images_compounded = num_images_compounded + 1;
        end
        image.data = compound_image / num_images_compounded;
        
        % envelope detect after if coherent compounding
        if strcmp(method,'coherent')
            image.data = envelopeDetection(image.data);
        end
    
	% recon without frequency compounding
    else
        
        % recon with gaussian bandpass filter
        if simu.params.gaussian_freq_filtered
            image.data = reconstruct2dUSimage(sensor.data, sensor.params, simu.params.c0, ...
                                        'Upsample', true, ...
                                        'FreqBandFilter', {simu.params.freq_filter_cf, simu.params.freq_filter_bw}, ...
                                        'EnvelopeDetect', true, ...
                                        'SaveImageToFile', false );
        
        % recon with custom frequency filter
        elseif simu.params.custom_freq_filtered
            image.data = reconstruct2dUSimage(sensor.data, sensor.params, simu.params.c0, ...
                                        'Upsample', true, ...
                                        'FreqCustomFilter', simu.params.custom_freq_filter_data, ...
                                        'EnvelopeDetect', true, ...
                                        'SaveImageToFile', false );
        
        % recon without frequency filter
        else
            image.data = reconstruct2dUSimage(sensor.data, sensor.params, simu.params.c0, ...
                                        'Upsample', true, ...
                                        'EnvelopeDetect', true, ...
                                        'SaveImageToFile', false );
        end
        
    end
    
    % assign image.kgrid and image.t_array from global
    global kgrid t_array
    image.kgrid   = kgrid;
    image.t_array = t_array;
    
end

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
        title([simu.params.scatt_type   ' c ' num2str(simu.params.scatt_c)  ' rho ' num2str(simu.params.scatt_rho) ', ' ...
               simu.params.object_shape ' c ' num2str(simu.params.object_c) ' rho ' num2str(simu.params.object_rho) ])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar
    
end

function save_image_for_sliceViewer(image, sensor, simu, image_file_path) %#ok<DEFNU>
% saves:    image.data reshaped to 3D
% requires: image.data
%           image.kgrid   - for volume spacing
%           sensor.params - for dt
%           simu.params   - for c0, file name

    [Nx_image, Ny_image] = size(image.data);
    volume_data = reshape(image.data, Nx_image, 1, Ny_image);
    volume_spacing = [image.kgrid.dx, 1, sensor.params.dt*simu.params.c0];
                     % omit factor 1/2 in dz because of doubled depth bug
    
    save([image_file_path '_image_4sliceViewer.mat'], 'volume_data', 'volume_spacing', '-v7.3')

end


%% IMAGE QUAL FUNCTIONS

function ROI = get_discROI_at_object_in_image(image, simu, fractional_radius)

    object_x = simu.kgrid.x_vec(simu.params.object_x);
    object_y = (simu.params.object_y - simu.params.pml_size - 1) * simu.kgrid.dx;
    object_radius = 0.5e-3;
    
    vec_x = image.kgrid.x_vec;
    vec_y = image.t_array*image.c0;
    
    distance = sqrt((vec_x - object_x).^2 + (vec_y - object_y).^2);
    
    ROI = distance < object_radius * fractional_radius;
    
end

function [scatter_hole_mean, scatter_hole_std] = get_scattering_distr_in_hole(image, simu, plot_toggle)

    mask = get_discROI_at_object_in_image(image, simu, 0.9);
    
    ROI = image.data(mask);
    
    scatter_hole_mean = mean(ROI(:));
    scatter_hole_std  = std(ROI(:));
    
    if plot_toggle == true
        plot_histogram_of_scattering_distr(ROI, 'scatter hole', 'Normalise', true)
    end
    
end

function [scatter_stmm_mean, scatter_stmm_std] = get_scattering_distr_in_stmm(image, simu, plot_toggle)

    hole_notoutside = get_discROI_at_object_in_image(image, simu, 1.1);
    hole_outside    = not(hole_notoutside);
    large_hole      = get_discROI_at_object_in_image(image, simu, 2);
    
    mask = and(hole_outside,large_hole);
    
    ROI = image.data(mask);
    
    scatter_stmm_mean = mean(ROI(:));
    scatter_stmm_std  = std(ROI(:));
    
    if plot_toggle == true
        plot_histogram_of_scattering_distr(ROI, 'scatter tmm', 'Normalise', true)
    end
    
end

function plot_histogram_of_scattering_distr(ROI, legend_entry, varargin)

    % set default
    num_req_input_variables = 2;
    toNormalise = false;
    
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'Normalise'
                    toNormalise = varargin{input_index + 1};
                otherwise
                    error('Unknown optional input.');
            end
        end
    end
    
    binwidth   = max(ROI(:))/100;    
    binedges   = 0:binwidth:max(ROI(:))*1.1+binwidth;
    bincount   = histcounts(ROI,binedges);
    bincentres = 0.5*(binedges(2:end)+binedges(1:end-1));
	
    if toNormalise
        bincount = bincount / max(bincount);
    end
    
    figure(gcf)
    % histogram(ROI,'BinWidth',10,'DisplayStyle','stairs','DisplayName',legend_entry)
    plot(bincentres,bincount,'DisplayName',legend_entry)
    
end

function signal_wire = get_specular_signal_of_wire(image, simu)

    mask = get_discROI_at_object_in_image(image, simu, 1);
    
    ROI = image.data(mask);
    
    signal_wire = max(ROI,[],'all');

end

function [resoLat, resoAxi] = get_resolution_of_wire(image, simu)

    signal_wire = get_specular_signal_of_wire(image, simu);
    
    [peak_pos_x, peak_pos_z] = find(image.data == signal_wire);
    
    halfwidth = 0.5e-3;
    hwX = round(halfwidth/image.kgrid.dx);
    hwZ = round(halfwidth/(simu.kgrid.dt*image.c0));
    
    profile_x = image.data(peak_pos_x-hwX:peak_pos_x+hwX, peak_pos_z);
    profile_z = image.data(peak_pos_x, peak_pos_z-hwZ:peak_pos_z+hwZ);
    
    resoLat = fwhm(profile_x,image.kgrid.dx,1);
    resoAxi = fwhm(profile_z,simu.kgrid.dt*image.c0,1);

end




