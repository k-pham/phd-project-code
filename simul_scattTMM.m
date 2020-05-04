% Pulse-echo plane-wave US imaging - scattering TMM with holes (simulate agarTMM)
% combine _rand and _points scripts into one


%% loops for parameter search

% for c_range = [0,1,2,4,10,30,50,100,150,300]
% for rho_range = [0,1,2,4,10,30,50,100,200]
%% SET UP

% clear all
close all

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';


%% SET UP SIMULATION -> struct SIMU

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% choose option to represent scattering medium as
simu.params.scattering_type = 'random';     % options: 'random', 'points'
simu.params.object_shape    = 'hole';       % options: 'hole', 'slab'

simu = set_simu_kgrid(simu);
simu = set_medium(simu, c0, rho0);
simu = make_time(simu, c0);
simu = set_pml(simu);
simu = set_source(simu);
simu = set_sensor(simu);
simu = set_inputs(simu);

fig_simu = plot_simu_medium(simu);


%% RUN SIMULATION -> struct SENSOR

sensor.data = kspaceFirstOrder2DC(simu.kgrid, simu.medium, simu.source, simu.sensor, simu.inputs{:});

sensor = set_params(sensor, simu);
sensor = set_sensor_kgrid(sensor, simu);
sensor = set_sensor_tarray(sensor, simu);

fig_sens = plot_sensor_data(sensor, simu);


%% RECONSTRUCTION -> struct IMAGE

image.data = reconstruct2dUSimage(sensor.data, sensor.params, c0);
    % NOTE: kgrid and t_array UPDATED
    global kgrid t_array

image.kgrid   = kgrid;
image.t_array = t_array;

fig_imag = plot_image_data(image, simu, c0);

save_image_for_sliceViewer(image, sensor, simu, file_dir_data, c0)


%% SAVE FIGURES & SIMU DATA & SENSOR DATA & IMAGE DATA

file_name = get_file_name(simu);

% saveas(fig_simu,[file_dir_figs file_name '_medium.fig'])
saveas(fig_simu,[file_dir_figs file_name '_medium.jpg'])
% saveas(fig_sens,  [file_dir_figs file_name '_data.fig'])
saveas(fig_sens,  [file_dir_figs file_name '_data.jpg'])
% saveas(fig_imag, [file_dir_figs file_name '_image.fig'])
saveas(fig_imag, [file_dir_figs file_name '_image.jpg'])

save([file_dir_data file_name '_simu']  , 'simu'  , '-v7.3')
save([file_dir_data file_name '_sensor'], 'sensor', '-v7.3')
save([file_dir_data file_name '_image'] , 'image' , '-v7.3')


%% end of loops for parameter search
% end
% end

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

function file_name = get_file_name(simu)

    file_name = [simu.params.scattering_type '_SCATT_c'  num2str(simu.params.c_scatt)  '_rho' num2str(simu.params.rho_scatt) '_' ...
                 simu.params.object_shape    '_OBJECT_c' num2str(simu.params.c_object) '_rho' num2str(simu.params.rho_object) ];

end


%% METHODS FOR SIMU

function simu = set_simu_kgrid(simu)

    dx = 10e-6;                 % grid point spacing in the x direction [m]
    dy = dx;                    % grid point spacing in the y direction [m]
    Nx = 1536;                  % number of grid points in the x (row) direction
    Ny = 1024;                  % number of grid points in the y (column) direction
    
    simu.kgrid = kWaveGrid(Nx, dx, Ny, dy);

end

function simu = set_medium(simu, c0, rho0)

    % define random scattering medium
    c_range   = 40;         % [m/s]
    rho_range = 80;         % [m/s]
    
    % define scattering medium with point scatterers
    c_pointscatt   = 1550;          % [m/s]
    rho_pointscatt = 1100;          % [m/s]
        num_points_per_voxel = 2;
        vox_size = 100e-6;          % [m]
    
    % define non-scattering object (hole/slab)
    c_object   = 1500;        % [m/s]
    rho_object = 1000;        % [m/s]
    
    % make background medium
    simu.medium.sound_speed = c0   * ones(simu.kgrid.Nx, simu.kgrid.Ny);
    simu.medium.density     = rho0 * ones(simu.kgrid.Nx, simu.kgrid.Ny);
    
    % make scattering medium specified by scattering_type
    switch simu.params.scattering_type
        case 'random'
            c_rand    = ( rand(simu.kgrid.Nx,simu.kgrid.Ny) - 0.5 ) * c_range;
            rho_rand  = ( rand(simu.kgrid.Nx,simu.kgrid.Ny) - 0.5 ) * rho_range;
            
            simu.medium.sound_speed = simu.medium.sound_speed + c_rand;
            simu.medium.density     = simu.medium.density     + rho_rand;
            
            simu.params.c_scatt   = c_range;
            simu.params.rho_scatt = rho_range;
        case 'points'
            pointscatts = get_pointscatt_locations(simu.kgrid.Nx, simu.kgrid.Ny, ...
                                                   simu.kgrid.dx, simu.kgrid.dy, ...
                                                   num_points_per_voxel, vox_size);
            
            simu.medium.sound_speed(pointscatts==1) = c_pointscatt;
            simu.medium.density(pointscatts==1)     = rho_pointscatt;
            
            simu.params.c_scatt   = c_pointscatt;
            simu.params.rho_scatt = rho_pointscatt;
    end
    
    % make object specified by object_shape
    switch simu.params.object_shape
        case 'hole'
            object = get_hole_location(simu.kgrid.Nx, simu.kgrid.Ny);
        case 'slab'
            object = get_slab_location(simu.kgrid.Nx, simu.kgrid.Ny);
    end
    
    simu.medium.sound_speed(object==1) = c_object;
    simu.medium.density(object==1)     = rho_object;
	
    simu.params.c_object   = c_object;
    simu.params.rho_object = rho_object;

end

function simu = make_time(simu, c0)

    cfl = 0.2;                  % clf number
    shorten_time = 0.75;        % fraction by which to shorten acquisition length
    t_end = shorten_time * 2 * simu.kgrid.Ny * simu.kgrid.dy / c0;      % [s]
    
    simu.kgrid.makeTime(simu.medium.sound_speed, cfl, t_end);

end

function simu = set_pml(simu)

    simu.params.pml_size = 20;              % thickness of the PML [grid points]

end

function simu = set_source(simu)

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

function simu = set_sensor(simu)

    sensor_spacing      = 100e-6;                                   % [m]
    sensor_spacing_grid = round(sensor_spacing / simu.kgrid.dx);    % [grid points]
    sensor_positions_x  = simu.params.pml_size : sensor_spacing_grid : (simu.kgrid.Nx - simu.params.pml_size);
    
    simu.sensor.mask = zeros(simu.kgrid.Nx, simu.kgrid.Ny);
    simu.sensor.mask(sensor_positions_x, simu.params.pml_size+1) = 1;
    
    simu.params.sensor_spacing_grid = sensor_spacing_grid;

end

function simu = set_inputs(simu)

    simu.inputs = {'PMLSize', simu.params.pml_size, 'PlotLayout', true, 'PlotSim', true};

end

function fig_simu = plot_simu_medium(simu)

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

function sensor = set_params(sensor, simu)

    sensor.params.Nx = size(sensor.data,1);
    sensor.params.Ny = 1;
    sensor.params.dx = simu.params.sensor_spacing_grid * simu.kgrid.dx;
    sensor.params.dy = simu.kgrid.dy;
    sensor.params.dt = simu.kgrid.dt;
    
    sensor.params.trigger_delay      = 0;
    sensor.params.Nt_zero_pad_source = 50;
    sensor.params.Nt_t0_correct      = -16;
    sensor.params.file_data          = '111111\scattTMM_simul';

end

function sensor = set_sensor_kgrid(sensor, simu)

    sensor.kgrid    = kWaveGrid(sensor.params.Nx, sensor.params.dx, sensor.params.Ny, sensor.params.dy);
    sensor.kgrid.dt = simu.kgrid.dt;
    sensor.kgrid.Nt = simu.kgrid.Nt;

end

function sensor = set_sensor_tarray(sensor, simu)

    sensor.t_array  = simu.kgrid.t_array;

end

function fig_sens = plot_sensor_data(sensor, simu)

    % figure
    % imagesc(sensor_data(:,1:50)')
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

function fig_imag = plot_image_data(image, simu, c0)

    fig_imag = figure;
    imagesc(image.kgrid.x_vec*1e3,image.t_array*c0*1e3,image.data')     % omit factor 1/2 in depth because of doubled depth bug
        axis image
        title([simu.params.scattering_type ' c ' num2str(simu.params.c_scatt) ' rho ' num2str(simu.params.rho_scatt)])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar

end

function save_image_for_sliceViewer(image, sensor, simu, file_dir_data, c0)

    [Nx_image, Ny_image] = size(image.data);
    volume_data = reshape(image.data, Nx_image, 1, Ny_image);
    volume_spacing = [image.kgrid.dx, 1, sensor.params.dt*c0];	% omit factor 1/2 in dz because of doubled depth bug
    
    file_name = get_file_name(simu);
    
    save([file_dir_data file_name '_image_4sliceViewer.mat'], 'volume_data', 'volume_spacing', '-v7.3')

end

