% Pulse-echo plane-wave US imaging - scattering TMM with holes (simulate agarTMM)
% combine _rand and _points scripts into one

global kgrid t_array

%% loops for parameter search

% for c_range = [0,1,2,4,10,30,50,100,150,300]
% for rho_range = [0,1,2,4,10,30,50,100,200]
%% SET UP

% clear all
close all

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';


%% SET UP SIMULATION

% create the computational grid
dx = 10e-6;                 % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
Nx = 1536;                  % number of grid points in the x (row) direction
Ny = 1024;                  % number of grid points in the y (column) direction
pml_size = 20;              % thickness of the PML [grid points]

% background material properties (WATER)
c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]

% choose option to represent scattering medium as
simu.medium.scattering_type = 'random';     % options: 'random', 'points'
simu.medium.object_shape    = 'hole';       % options: 'hole', 'slab'

simu.kgrid = kWaveGrid(Nx, dx, Ny, dy);
simu = set_medium(simu, c0, rho0);
simu = make_time(simu, c0, 0.75);
simu = set_source(simu, pml_size);
simu = set_sensor(simu, pml_size);
simu.inputs = {'PMLSize', pml_size, 'PlotLayout', true, 'PlotSim', true};

fig_simu = plot_simu_medium(simu);


%% RUN SIMULATION

sensor_data = kspaceFirstOrder2DC(simu.kgrid, simu.medium, simu.source, simu.sensor, simu.inputs{:});

% figure
% imagesc(sensor_data(:,1:50)')
%     xlabel('x position / dx')
%     ylabel('time / dt')

fig_data = figure;
imagesc(simu.kgrid.x_vec*1e3,simu.kgrid.t_array(50:end)*1e6,sensor_data(:,50:end)')
    title([scattering_type ' c ' num2str(simu.medium.c_scatt) ' rho ' num2str(simu.medium.rho_scatt)])
    xlabel('x position / mm')
    ylabel('time / \mus')
    colorbar

params.Nx = size(sensor_data,1);
params.Ny = 1;
params.dx = 10*dx;
params.dy = dy;
params.dt = simu.kgrid.dt;

params.trigger_delay      = 0;
params.Nt_zero_pad_source = 50;
params.Nt_t0_correct      = -16;
params.file_data          = '111111\scattTMM_simul';


%% IMAGE FORMATION

reflection_image = reconstruct2dUSimage(sensor_data, params, c0);
    % NOTE: kgrid and t_array UPDATED

fig_image = figure;
% imagesc(kgrid.x_vec*1e3,kgrid.t_array*c0/2*1e3,reflection_image')
%     % can use kgrid.x_vec here, even though spacing is different in image,
%     % because only care about end points of sensor; kgrid.t_array is
%     % correct to use with c0/2 scaling
imagesc(kgrid.x_vec*1e3,t_array*c0*1e3,reflection_image')
    % using updated kgrid and t_array from reconstruct2dUSimage
    axis image
    title([scattering_type ' c ' num2str(simu.medium.c_scatt) ' rho ' num2str(simu.medium.rho_scatt)])
    xlabel('x position / mm')
    ylabel('y position / mm')
    colorbar

% prepare for saving
[Nx_image, Ny_image] = size(reflection_image);
volume_data = reshape(reflection_image,Nx_image,1,Ny_image);
volume_spacing = [kgrid.dx, 1, params.dt*c0];           % omit factor 1/2 in dz because of doubled depth bug


%% SAVE FIGURES & SENSOR DATA & IMAGE DATA

% file_name = [scattering_type '_SCATT_c' num2str(simu.medium.c_scatt) '_rho' num2str(simu.medium.rho_scatt) ...
%              object_shape    '_OBJECT_c' num2str(simu.medium.c_object) '_rho' num2str(simu.medium.rho_object) ];
% 
% % saveas(fig_medium,[file_dir_figs file_name '_medium.fig'])
% saveas(fig_medium,[file_dir_figs file_name '_medium.jpg'])
% % saveas(fig_data,  [file_dir_figs file_name '_data.fig'])
% saveas(fig_data,  [file_dir_figs file_name '_data.jpg'])
% % saveas(fig_image, [file_dir_figs file_name '_image.fig'])
% saveas(fig_image, [file_dir_figs file_name '_image.jpg'])
% 
% save([file_dir_data file_name '_sensor_data'], 'sensor_data', 'params', 'simu')
% 
% save([file_dir_data file_name '_image_data.mat'], 'volume_data', 'volume_spacing', 'kgrid', 't_array', '-v7.3')


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


%% METHODS FOR SIMU

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
    switch simu.medium.scattering_type
        case 'random'
            c_rand    = ( rand(simu.kgrid.Nx,simu.kgrid.Ny) - 0.5 ) * c_range;
            rho_rand  = ( rand(simu.kgrid.Nx,simu.kgrid.Ny) - 0.5 ) * rho_range;
            
            simu.medium.sound_speed = simu.medium.sound_speed + c_rand;
            simu.medium.density     = simu.medium.density     + rho_rand;
            
            simu.medium.c_scatt   = c_range;
            simu.medium.rho_scatt = rho_range;
        case 'points'
            pointscatts = get_pointscatt_locations(simu.kgrid.Nx, simu.kgrid.Ny, ...
                                                   simu.kgrid.dx, simu.kgrid.dy, ...
                                                   num_points_per_voxel, vox_size);
            
            simu.medium.sound_speed(pointscatts==1) = c_pointscatt;
            simu.medium.density(pointscatts==1)     = rho_pointscatt;
            
            simu.medium.c_scatt   = c_pointscatt;
            simu.medium.rho_scatt = rho_pointscatt;
    end
    
    % make object specified by object_shape
    switch simu.medium.object_shape
        case 'hole'
            object = get_hole_location(simu.kgrid.Nx, simu.kgrid.Ny);
        case 'slab'
            object = get_slab_location(simu.kgrid.Nx, simu.kgrid.Ny);
    end
    
    simu.medium.sound_speed(object==1) = c_object;
    simu.medium.density(object==1)     = rho_object;
	
    simu.medium.c_object   = c_object;
    simu.medium.rho_object = rho_object;

end

function simu = make_time(simu, c0, shorten_time)

    cfl = 0.2;
    t_end = shorten_time*2*simu.kgrid.Ny*simu.kgrid.dy/c0;      % [s]
    
    simu.kgrid.makeTime(simu.medium.sound_speed,cfl,t_end);

end

function simu = set_source(simu, pml_size)

    amplitude = 10;             % [Pa]
    
    apodisation_width = 0.4;    % proportion of width Nx
    apodisation = getWin(simu.kgrid.Nx, 'Gaussian', 'Param', apodisation_width);
    
    pulse_tpeak      = 20e-9;   % time of pulse peak pressure [s]
    pulse_width      = 16e-9;   % FWHM-width of pulse [s]
    pulse_variance   = (pulse_width / ( 2*sqrt(2*log(2)) ) )^2;
    pulse = gaussian(simu.kgrid.t_array, 1, pulse_tpeak, pulse_variance);
    
    simu.source.p_mask = zeros(simu.kgrid.Nx, simu.kgrid.Ny);
    simu.source.p_mask(:, pml_size + 1) = 1;
    
    simu.source.p = amplitude * apodisation * pulse;

end

function simu = set_sensor(simu, pml_size)

    sensor_spacing      = 100e-6;   % [m]
    sensor_spacing_grid = round(sensor_spacing / simu.kgrid.dx);
    sensor_positions_x  = pml_size : sensor_spacing_grid : (simu.kgrid.Nx - pml_size);
    
    simu.sensor.mask = zeros(simu.kgrid.Nx, simu.kgrid.Ny);
    simu.sensor.mask(sensor_positions_x, pml_size+1) = 1;

end

function fig_simu = plot_simu_medium(simu)

    fig_simu = figure;
    set(gcf,'Position',[200,200,500,700])
    subplot(2,1,1)
    imagesc(simu.kgrid.x_vec*1e3,simu.kgrid.y_vec*1e3,simu.medium.sound_speed')
        axis image
        title(['sound speed: ' simu.medium.scattering_type num2str(simu.medium.c_scatt)])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar
    subplot(2,1,2)
    imagesc(simu.kgrid.x_vec*1e3,simu.kgrid.y_vec*1e3,simu.medium.density')
        axis image
        title(['density: ' simu.medium.scattering_type num2str(simu.medium.rho_scatt)])
        xlabel('x position / mm')
        ylabel('y position / mm')
        colorbar

end

