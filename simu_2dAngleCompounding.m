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
medium.sound_speed = medium.sound_speed + c_rand;
medium.density     = medium.density     + rho_rand;
medium.sound_speed(hole==1) = hole_c;
medium.density(hole==1)     = hole_rho;

% create the time array
cfl = 0.2;                  % CFL number
shorten_time = 0.6;         % fraction to shorten length of time array
t_end = shorten_time*2*Ny*dy/c0;     % end time of the simulation [s]
kgrid.makeTime(medium.sound_speed,cfl,t_end);

% -------------------------------------------------------------------------
%                   make *angled plane wave* source
% -------------------------------------------------------------------------
%% anle loop
for source_angle = [-20:5:20,-14:2:14]
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
% sensor_data_sourceOnly = sensor_data;
% save(file_data_sourceOnly,'sensor_data_sourceOnly')
% end


%% ========================================================================
%                        PROCESS DATA (T_0 & ANGLE)
% =========================================================================

%% extract TOA from timeseries data

[~,TOA_x_idx] = max(sensor_data,[],2);     % [time index]
TOA_x = TOA_x_idx * kgrid.dt;              % [s]


%% fit plane wave (line) to source

% restrict fit range to exclude outliers at the edge
fit_mask = 62:282;

% regression in 1d TOA(x) = p1 + p2*x
Model = [ones(length(sensor_kgrid.x_vec(fit_mask)),1) sensor_kgrid.x_vec(fit_mask)];
params = Model \ TOA_x(fit_mask);

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


%% subtract acoustic source from data using saved data

sensor_data_sourceOnly = load(file_data_sourceOnly,'sensor_data_sourceOnly');
sensor_data_sourceOnly = sensor_data_sourceOnly.sensor_data_sourceOnly;

sensor_data = sensor_data - sensor_data_sourceOnly;

% plot the source-subtracted data
figure
imagesc(kgrid.t_array*1e6,sensor_positions,sensor_data)
xlabel('time [\mus]')
ylabel('sensor number')


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
for t0_excitation = floor(t0_excitation_true)%3%-500:200:500 % [dt]


%% zero padding source offset to sensor plane
% assumes that correction involves zero padding rather than trimming

% set t0(source) in dt
t0_source = t0_excitation - (TOA_source - t0_excitation);

% number of pads
pads = -t0_source;
% pads = -t0_excitation;

if pads >= 0
    sensor_data_padded = [ zeros(sensor_kgrid.Nx,pads) sensor_data ];
elseif pads < 0
    sensor_data_padded = sensor_data(:,-pads+1:end);
end


%% reconstruction

reflection_image = kspaceLineRecon_US_steered(sensor_data_padded',sensor_kgrid.dx,kgrid.dt,c0,a);
% reflection_image = kspaceLineRecon_US(sensor_data_padded',sensor_kgrid.dx,kgrid.dt,c0);
                                                        % input p_tx, output p_zx
reflection_image = permute(reflection_image,[2 1]);     % reorder p_zx to p_xz
reflection_image = reflection_image(:,1:round(kgrid.Nt/2));  % trim half of z if using doubled depth kgrid


%% envelope detection

disp('Envelope detecting ...')
tic
% reflection_image_env = envelopeDetection(reflection_image);
reflection_image_env = transpose(envelopeDetection(reflection_image'));
assert(isequal( size(reflection_image_env), size(reflection_image) ))
disp(['  completed in ' scaleTime(toc)]);


%% plot image

y_vec = kgrid.t_array * c0;
y_vec = y_vec(1:round(kgrid.Nt/2));

fig_img = figure;
%imagesc(reflection_image_env')
imagesc(sensor_kgrid.x_vec*1e3,y_vec*1e3,reflection_image_env')
title(['angle ' num2str(source_angle) ', t0 = ' num2str(t0_excitation) '*dt'])
ylabel('depth / mm')
xlabel('x axis / mm')
colorbar
axis image


%% save sensor & image data and images

file_specifier = ['_' num2str(dx*1e6) 'um' ...
                  '_' num2str(shorten_time) 't_end' ...
                  '_offsetSource' num2str(source_offset_y*1e3) 'e-3' ...
                  '_angle' num2str(source_angle) ...
                  '_scatt_c' num2str(scatt_c) '_rho' num2str(scatt_rho) ];

file_data_image = [file_dir_data 'scattTMM' file_specifier];
file_figs_image = [file_dir_figs 'scattTMM' file_specifier];

save([file_data_image '_data.mat'],'kgrid','sensor_kgrid','sensor_data','sensor_data_padded','y_vec','reflection_image','reflection_image_env');
% save([file_data_image '_image_data_raw.mat']      ,'reflection_image'    );
% save([file_data_image '_image_data_envDetect.mat'],'reflection_image_env');
% save([file_data_image '_image_data_params.mat']   ,'kgrid','sensor_kgrid','y_vec');

savefig(fig_img,[file_figs_image '.fig'])
saveas( fig_img,[file_figs_image '.jpg'])


%%
% pause
end % t0


%% ========================================================================
%                              IMAGE QUALITY
% =========================================================================

%% get hole location in image

hole_x_img = kgrid.x_vec(hole_x);                   % x position [m]
hole_y_img = (hole_y - pml_size -1) * kgrid.dy;     % y position [m]
hole_r_img = hole_radius * kgrid.dx;                % radius [m]

distance = sqrt((sensor_kgrid.x_vec - hole_x_img).^2 + (y_vec - hole_y_img).^2);    % [m]

mask_inside     = distance < 0.8 * hole_r_img;
mask_notoutside = distance < 1.2 * hole_r_img;
mask_outside    = not(mask_notoutside);
mask_largehole  = distance < 2   * hole_r_img;
mask_outerring  = and(mask_outside,mask_largehole);

% figure,imagesc(sensor_kgrid.x_vec*1e3,y_vec*1e3,mask_inside')   , axis image
% figure,imagesc(sensor_kgrid.x_vec*1e3,y_vec*1e3,mask_outerring'), axis image

ROI_hole = reflection_image_env(mask_inside);
ROI_stmm = reflection_image_env(mask_outerring);

scatter_hole_mean = mean(ROI_hole(:));
scatter_hole_std  = std(ROI_hole(:));

scatter_stmm_mean = mean(ROI_stmm(:));
scatter_stmm_std  = std(ROI_stmm(:));

scatSNR = scatter_stmm_mean / scatter_hole_mean;
scatCNR = (scatter_stmm_mean - scatter_hole_mean) / (scatter_stmm_std + scatter_hole_std);

figure, hold on
plot_histogram_of_scattering_distr(ROI_hole, 'hole', 'Normalise', true)
plot_histogram_of_scattering_distr(ROI_stmm, 'stmm', 'Normalise', true)
title(['scattering distributions, SNR = ' num2str(scatSNR) ', CNR = ' num2str(scatCNR)])
xlabel('pixel intensity')
ylabel('count')
legend(gca,'show')

disp('  scatSNR   scatCNR')
disp([scatSNR,scatCNR])

end


%% FUNCTIONS

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
