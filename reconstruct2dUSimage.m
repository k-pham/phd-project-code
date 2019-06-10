function [reflection_image, t_array, kgrid] = reconstruct2dUSimage(sensor_data, params, c0)


%% assign parameters from params to use in script and globally

global Nx Ny Nt dx dy dt

% create the computational grid
Nx = params.Nx;                         % number of grid points in the x (row) direction
Ny = params.Ny;                         % number of grid points in the y (column) direction - 1 for LINE SCAN
dx = params.dx;                         % grid point spacing in the x direction [m]
dy = params.dy;                         % grid point spacing in the y direction [m]

% scope parameters to make time array
Nt = size(sensor_data,2);               % number of samples in acquisition - alternatively: params.Nt
dt = params.dt;                         % dt between samples [s]
Nt_delay = params.trigger_delay / dt;   % number of samples to delay OR READ OUT FROM FILE NAME - HOW?

% if line scan along y, redefine as x axis
if Ny > Nx
    Nx = Ny;
    dx = dy;
end

%%

Nt = Nt_delay + Nt + params.Nt_t0_correct;
t_array = linspace(1,Nt,Nt)*dt;


%% prepare data for reconstruction

% remove the source from the time series & add zero padding for delay
sensor_data = [ zeros(Nx,int32(Nt_delay)) zeros(Nx,params.Nt_zero_pad_source) sensor_data(:,params.Nt_zero_pad_source+1:end) ];

% add/remove samples from sensor_data for t0 correction
if params.Nt_t0_correct > 0
    sensor_data = [ zeros(Nx,int32(params.Nt_t0_correct)) sensor_data ];
elseif params.Nt_t0_correct < 0
    sensor_data = sensor_data(:,-params.Nt_t0_correct+1:end);
end

%% upsample along space (x)
sensor_data_upsampled = zeros( 2*Nx , Nt );
for i = 1:Nx
    sensor_data_upsampled(2*i-1,:) = sensor_data(i,:);
    sensor_data_upsampled(2*i  ,:) = sensor_data(i,:);
end
sensor_data = sensor_data_upsampled;
clear sensor_data_upsampled
dx = dx/2;
Nx = 2*Nx;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% window the data (apodising to remove artefacts due to edge of sensor)
% win = getWin(Nx, 'Cosine');
% sensor_data_apodised = bsxfun(@times, win, sensor_data);

%% reconstruct an image using a k-space method
reflection_image = kspaceLineRecon_US(sensor_data_apodised', dx, dt, c0);   % input p_tx, output p_zx
reflection_image = permute(reflection_image,[2 1]);                         % reorder p_zx to p_xz

%% time gain compensation
% tgc_exponent = 50;
% tgc = 2+ exp(tgc_exponent * t_array * c0); % exponential
% %tgc = tgc_exponent * t_array * c0; % linear
% reflection_image = bsxfun(@times, tgc, reflection_image);

%% envelope detection
reflection_image = envelopeDetection(squeeze(reflection_image));            % p_xz

%% log compression
% compression_ratio = 3;
% reflection_image = logCompression(reflection_image, compression_ratio, true);   % p_xz

%% saving image data to .mat - NOT TESTED
[Nz] = size(reflection_image,2);
volume_data = reshape(reflection_image,Nx,1,Nz);
volume_spacing = [kgrid.dx, kgrid.dx, params.dt*c0];

phantom_id = strtok(params.file_data,'@'); % parse string up to specified delimiter
phantom_id = phantom_id(8:end);     % remove date folder from string
save(['recon_data\' phantom_id '.mat'],'volume_data','volume_spacing','-v7.3')



end
