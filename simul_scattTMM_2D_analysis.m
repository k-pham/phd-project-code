%% load images

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';

scattering_type = 'random';      % options: 'random', 'points'

c_ranges   = 0:10:150;
rho_ranges = 0:10:100;

c_hole     = 1500;
rho_hole   = 1000;

for idx_c = 16:16%length(c_scatts)
    for idx_r = 11:11%length(rho_scatts)

        c_scatt = c_ranges(idx_c);
        rho_scatt = rho_ranges(idx_r);
        
        file_name = [scattering_type '_SCATT_c' num2str(c_scatt) '_rho' num2str(rho_scatt) ...
                                     '_HOLE_c' num2str(c_hole) '_rho' num2str(rho_hole) ];
        
        SNS = load([file_dir_data file_name '_sensor_data.mat']);
        IMG = load([file_dir_data file_name '_image_data.mat']);
        
        sensor.data   = SNS.sensor_data;
        sensor.params = SNS.params;
        simu          = SNS.simu;
        
        sensor.kgrid   = kWaveGrid(sensor.params.Nx, sensor.params.dx, sensor.params.Ny, sensor.params.dy);
        sensor.t_array = simu.kgrid.t_array;
        
        image.data    = squeeze(IMG.volume_data);
        image.spacing = IMG.volume_spacing;
        image.kgrid   = IMG.kgrid;
        image.t_array = IMG.t_array;
        
        look_at_central_sensor_data(sensor)

    end
end


%% LOCAL FUNCTIONS

function look_at_central_sensor_data(sensor)

    x_centre = round(sensor.params.Nx/2);
    
    figure
    plot(sensor.t_array*1e6,sensor.data(x_centre,:)') % (50:end)
    
    figure
    imagesc(sensor.kgrid.x_vec*1e3,sensor.t_array(50:end)*1e6,sensor.data(:,50:end)')    

end

function image_quality_metric = get_image_quality_metric(reflection_image)

    

end





