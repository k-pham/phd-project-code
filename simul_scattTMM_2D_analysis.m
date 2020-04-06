%% load images

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';

scattering_type = 'random';      % options: 'random', 'points'

c_ranges   = 0:10:150;
rho_ranges = 0:10:100;

c0 = 1500;
rho0 = 1000;
c_hole     = 1500;
rho_hole   = 1000;

for idx_c = 5%length(c_scatts)
    for idx_r = 9%length(rho_scatts)

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
        
        % look_at_central_sensor_data(sensor)
        
        hole             = get_hole_location_in_image(image, c0, 1);
        hole_inside      = get_hole_location_in_image(image, c0, 0.9);
        hole_notoutside  = get_hole_location_in_image(image, c0, 1.1);
        
        hole_outside     = not(hole_notoutside);
        large_hole       = get_hole_location_in_image(image,c0, 2);
        
        hole_surrounding = and(hole_outside,large_hole);
        
        % figure
        % imagesc(image.kgrid.x_vec*1e3,image.t_array*c0*1e3,image.data')
        %     axis image
        %     title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt)])
        %     xlabel('x position / mm')
        %     ylabel('y position / mm')
        %     colorbar
        % 
        % figure
        % imagesc(image.kgrid.x_vec*1e3,image.t_array*c0*1e3,hole_surrounding')
        %     axis image
        %     title([scattering_type ' c ' num2str(c_scatt) ' rho ' num2str(rho_scatt)])
        %     xlabel('x position / mm')
        %     ylabel('y position / mm')
        %     colorbar
        
        

    end
end


%% compare with experimental data
% 181204 atmm with orgasol

file_dir = '../data/imagingUS/';
file_name = '181204/atmm_orgasol1_BK31[CNT]@0nm_t0[0]_dx[100�m]_dy[100�m]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
trigger_delay = 0;
samples_cut_off = 50;
samples_t0_correct = -6;
c0 = 1544;

[exp.data, exp.params] = loadSGL([file_dir file_name]);
exp.kgrid   = kWaveGrid(exp.params.Nx, exp.params.dx, exp.params.Ny ,exp.params.dy);
exp.t_array = linspace(1, exp.params.Nt, exp.params.Nt) * exp.params.dt;

figure
imagesc(exp.kgrid.x_vec*1e3, exp.t_array*1e6, squeeze(exp.data(:,round(exp.params.Nx/2),:))')
    xlabel('x / mm')
    ylabel('time / \mus')
    ylim([0,10])
figure
plot(exp.t_array*1e6, squeeze(exp.data(round(exp.params.Nx/2),round(exp.params.Ny/2),:)))
    xlabel('time / \mus')
    ylabel('amplitude')
    xlim([0,10])
    ylim([-0.3,0.2])


%% LOCAL FUNCTIONS

function look_at_central_sensor_data(sensor)

    x_centre = round(sensor.params.Nx/2);
    
    figure
    plot(sensor.t_array*1e6,sensor.data(x_centre,:)) % (50:end)
        xlabel('time / \mus')
        ylabel('acoustic pressure / Pa')
    
    % figure
    % imagesc(sensor.kgrid.x_vec*1e3,sensor.t_array(50:end)*1e6,sensor.data(:,50:end)')    

end

function hole = get_hole_location(Nx, Ny)

    hole_radius = 50;
    hole_x     = round(Nx/2);
    hole_y     = round(Ny/4);
    hole = makeDisc(Nx,Ny,hole_x,hole_y,hole_radius);
        
end

function hole = get_hole_location_in_image(image, c0, rad_frac)

    hole_x = 0;
    hole_y = 2.5e-3;
    hole_radius = 0.5e-3;
    
    vec_x = image.kgrid.x_vec;
    vec_y = image.t_array*c0;
    
    distance = sqrt((vec_x - hole_x).^2 + (vec_y - hole_y).^2);
    
    hole = distance < hole_radius * rad_frac;
    
end

function image_quality_metric = get_image_quality_metric(reflection_image)

    

end





