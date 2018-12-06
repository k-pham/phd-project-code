function [reflection_image, samples_total, t_array, kgrid] = ...
    reconstruct3dUSimage(sensor_data, params, trigger_delay, samples_cut_off, samples_t0_correct, c0, varargin)

% set usage defaults
num_req_input_variables = 6;
zero_pad_sides = 0;
toUpsample = false;
toApodise = false;
toTimeGainCompensate = false;
toEnvelopeDetect = false;
toLogCompress = false;


% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'ZeroPad'
                zero_pad_sides = varargin{input_index + 1};
            case 'Upsample'
                toUpsample = varargin{input_index + 1};
            case 'Apodise'
                toApodise = varargin{input_index + 1};
            case 'TimeGainCompensate'
                toTimeGainCompensate = varargin{input_index + 1};
            case 'EnvelopeDetect'
                toEnvelopeDetect = varargin{input_index + 1};
            case 'LogCompress'
                toLogCompress = varargin{input_index + 1};
            otherwise
                error('Unknown optional input.');
        end
    end
end

%% define parameters

% create the computational grid
Nx = params.Nx;                      % number of grid points in the x (row) direction
Ny = params.Ny;                      % number of grid points in the y (column) direction
dx = params.dx;                      % grid point spacing in the x direction [m]
dy = params.dy;                      % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% scope parameters to make time array
dt = params.dt;                         % dt between samples [s]
%trigger_delay = 8e-6;                  % delay of acquisition post-trigger [s]
samples_delay = trigger_delay / dt;     % number of samples to delay OR READ OUT FROM FILE NAME - HOW?
samples_data = size(sensor_data,3);     % number of samples in acquisition
%samples_t0_correct = 0;                % number of samples to add to correct for t0
samples_total = samples_delay + samples_data + samples_t0_correct;
t_array = linspace(1,samples_total,samples_total)*dt;

% medium parameters
%c0 = 1467;                          % sound speed [m/s]


%% image reconstruction

% remove the source from the time series
sensor_data = cat(3, zeros(Nx,Ny,samples_cut_off), sensor_data(:,:,samples_cut_off+1:end) );

% add zero padding for delay
sensor_data = cat(3, zeros(Nx,Ny,int32(samples_delay)), sensor_data );

% add/remove samples from sensor_data for t0 correction
if samples_t0_correct > 0
    sensor_data = cat(3, zeros(Nx,Ny,int32(samples_t0_correct)), sensor_data );
elseif samples_t0_correct < 0
    sensor_data = sensor_data(:,:,-samples_t0_correct+1:end);
end

% zero pad the sides
if zero_pad_sides
    sensor_data = zero_padding_sides(sensor_data,zero_pad_sides);
end

% upsample along space (x and y)
if toUpsample
    sensor_data_upsampled = zeros( 2*Nx , 2*Ny, samples_total );
    for i = 1:Nx
        for j = 1:Ny
            sensor_data_upsampled(2*i-1,2*j-1,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i-1,2*j  ,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i  ,2*j-1,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i  ,2*j  ,:) = sensor_data(i,j,:);
        end
    end
    sensor_data = sensor_data_upsampled;
    clear sensor_data_upsampled
    dx = dx/2;
    dy = dy/2;
    Nx = 2*Nx;
    Ny = 2*Ny;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
end

% window the data (apodising to remove artefacts due to edge of sensor)
% win = getWin([Nx Ny], 'Cosine');
% win = win + 0.5;
% sensor_data_apodised = bsxfun(@times, win, sensor_data);
sensor_data_apodised = sensor_data;

% reconstruct an image using a k-space method
sensor_data_apodised = permute(sensor_data_apodised,[3 1 2]);       % reorder p_xyt to p_txy
reflection_image = kspacePlaneRecon_US(sensor_data_apodised, dx, dy, dt, c0);   % output as p_zxy
reflection_image = permute(reflection_image,[2 3 1]);               % reorder p_zxy to p_xyz

% % time gain compensation
% tgc_exponent = 0.3;
% % tgc = exp(tgc_exponent * t_array * c0); % exponential
% tgc = tgc_exponent * t_array * c0; % linear
% tgc = reshape(tgc, 1, 1, length(tgc));
% reflection_image = bsxfun(@times, tgc, reflection_image);

% % positivity constraint (remove all negative intensities in image, e.g. back of sensor) - DOESN'T WORK SO WELL
% reflection_image(reflection_image<0) = 0;

% envelope detection per slice
for i = 1:Ny
    reflection_image(:,i,:) = envelopeDetection(squeeze(reflection_image(:,i,:)));
end

% % log compression
% compression_ratio = 3;
% reflection_image = logCompression(reflection_image, compression_ratio, true);

end


function sensor_data_padded = zero_padding_sides(sensor_data, pads)

    sensor_data_padded = zeros(Nx+2*pads, Ny+2*pads, samples_total);
    sensor_data_padded(pads+1:Nx+pads, pads+1:Ny+pads, :) = sensor_data;
    Nx = Nx + 2*pads;
    Ny = Ny + 2*pads;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    
end


%% plot slice of the reconstructed image
% 
% slice = 70;
% reflection_image_slice = squeeze(reflection_image(:,slice,1:samples_total/2)); % slice of p_xz
% 
% figure(7); clf(7)
% imagesc(t_array*c0/2*1e3, kgrid.x_vec*1e3, reflection_image_slice)
% % 1st index (x) = row index = y axis
%     title(['slice y = ' num2str(slice) ' through the reconstructed image'])
%     xlabel('z [mm]')
%     ylabel('x [mm]')
%     caxis([0,200])
%     %axis image
%     cmap = colormap(gray);
%     cmap = flipud(cmap);    % flip colormap to make black = signal
%     colormap(cmap);
%     colorbar
%     hold on
%     plot([0,0],[-5,5],'b--')
%     legend(['sensor'])
%     hold off
%     
% 
%% plot MIP of the reconstructed image
% 
% reflection_image_MIP = squeeze(max(reflection_image(:,:,:),[],3));      % p_xyz to p.max(z)_xy
% 
% figure(8); clf(8)
% imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, reflection_image_MIP')
% % need to transpose so that first index = y = row index = y-axis
%     title('MIP of reconstructed image')
%     xlabel('x [mm]')
%     ylabel('y [mm]')
%     axis image
%     colormap(gray)
%     cmap = colormap;
% %     cmap = flipud(cmap);
%     colormap(cmap);
%     colorbar
% 
%     
%% fly through reconstructed image
% 
% figure(9)
% for i = 1:1600
%     imagesc(kgrid.x_vec*1e3, kgrid.y_vec*1e3, squeeze(reflection_image(:,:,i))') % need to transpose ???
%         axis equal
%         colormap(gray)
%         cmap = colormap;
%         cmap = flipud(cmap);
%         colormap(cmap);
%         colorbar
%         caxis([-400 400])
%         title(i)
%         pause
%       
% end
%      
    
