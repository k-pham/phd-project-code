%% parameters needed to make kgrid and time

dx = 10e-6;                 % grid point spacing in the x direction [m]
dy = dx;                    % grid point spacing in the y direction [m]
Nx = 1536;                  % number of grid points in the x (row) direction
Ny = 1024;                  % number of grid points in the y (column) direction
kgrid = kWaveGrid(Nx, dx, Ny, dy);

c0 = 1500;      % sound speed [m/s]
rho0 = 1000;    % density [kg/m^3]


medium = define_random_medium(Nx,Ny,c0,rho0,c_range,rho_range,c_hole,rho_hole);


%% load images

file_dir_data = 'D:\PROJECT\data\simulations\scattTMM\';
file_dir_figs = 'D:\PROJECT\figures\_Matlab figs\simulations\scattTMM\';

scattering_type = 'random';      % options: 'random', 'points'

c_scatts   = 0:10:150;
rho_scatts = 0:10:100;

c_hole     = 1500;
rho_hole   = 1300;

for idx_c = 16:16%length(c_scatts)
    for idx_r = 11:11%length(rho_scatts)

        c_scatt = c_scatts(idx_c);
        rho_scatt = rho_scatts(idx_r);
        
        file_name = [scattering_type '_SCATT_c' num2str(c_scatt) '_rho' num2str(rho_scatt) ...
                                     '_HOLE_c' num2str(c_hole) '_rho' num2str(rho_hole) ];
        
        SNS = load([file_dir_data file_name '_sensor_data.mat']);
        IMG = load([file_dir_data file_name '_image_data.mat']);
        
        sensor_data = SNS.sensor_data;
        params      = SNS.params;
        
        image_data    = squeeze(IMG.volume_data);
        image_spacing = IMG.volume_spacing;
        

    end
end

