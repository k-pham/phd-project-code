function [reflection_image, samples_total, t_array, kgrid] = ...
    reconstruct2dUSimage(sensor_data, params, trigger_delay, samples_cut_off, samples_t0_correct, c0)


%% define parameters

% create the computational grid
Nx = params.Nx;                         % number of grid points in the x (row) direction
Ny = params.Ny;                         % number of grid points in the y (column) direction - 1 for LINE SCAN
dx = params.dx;                         % grid point spacing in the x direction [m]
dy = params.dy;                         % grid point spacing in the y direction [m]
%kgrid = kWaveGrid(Nx, dx, Ny, dy);

% scope parameters to make time array
dt = params.dt;                         % dt between samples [s]
%trigger_delay = 0;                     % delay of acquisition post-trigger [s]
samples_delay = trigger_delay / dt;     % number of samples to delay OR READ OUT FROM FILE NAME - HOW?
samples_data = size(sensor_data,2);     % number of samples in acquisition
%samples_t0_correct = 0;                % number of samples to add to correct for t0
samples_total = samples_delay + samples_data + samples_t0_correct;
t_array = linspace(1,samples_total,samples_total)*dt;

% medium parameters
%c0 = 1484;                             % sound speed [m/s]


%% image reconstruction

% remove the source from the time series & add zero padding for delay
sensor_data = [ zeros(Nx,int32(samples_delay)) zeros(Nx,samples_cut_off) sensor_data(:,samples_cut_off+1:end) ];

% add/remove samples from sensor_data for t0 correction
if samples_t0_correct > 0
    sensor_data = [ zeros(Nx,int32(samples_t0_correct)) sensor_data ];
elseif samples_t0_correct < 0
    sensor_data = sensor_data(:,-samples_t0_correct+1:end);
end

% upsample along space (x)
sensor_data_upsampled = zeros( 2*Nx , samples_total );
for i = 1:Nx
    sensor_data_upsampled(2*i-1,:) = sensor_data(i,:);
    sensor_data_upsampled(2*i  ,:) = sensor_data(i,:);
end
sensor_data = sensor_data_upsampled;
clear sensor_data_upsampled
dx = dx/2;
Nx = 2*Nx;
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% window the data (apodising to remove artefacts due to edge of sensor)
win = getWin(Nx, 'Cosine');
sensor_data_apodised = bsxfun(@times, win, sensor_data);

% reconstruct an image using a k-space method
reflection_image = kspaceLineRecon_US(sensor_data_apodised', dx, dt, c0);   % input p_tx, output p_zx
reflection_image = permute(reflection_image,[2 1]);                         % reorder p_zx to p_xz

% time gain compensation
% tgc_exponent = 50;
% tgc = 2+ exp(tgc_exponent * t_array * c0); % exponential
% %tgc = tgc_exponent * t_array * c0; % linear
% reflection_image = bsxfun(@times, tgc, reflection_image);

% envelope detection
reflection_image = envelopeDetection(squeeze(reflection_image));            % p_xz

% log compression
% compression_ratio = 3;
% reflection_image = logCompression(reflection_image, compression_ratio, true);   % p_xz


end
