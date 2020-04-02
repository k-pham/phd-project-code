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
        
        sensor_data   = SNS.sensor_data;
        sensor_params = SNS.params;
        simu          = SNS.simu;
        
        image_data    = squeeze(IMG.volume_data);
        image_spacing = IMG.volume_spacing;
        image_kgrid   = IMG.kgrid;
        image_t_array = IMG.t_array;
        

    end
end


%% LOCAL FUNCTIONS

function look_at_central_sensor_data(sensor_data, params)

    

end

function image_quality_metric = get_image_quality_metric(reflection_image)

    

end





