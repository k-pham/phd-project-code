function [reflection_image, t_array, kgrid] = ...
    reconstruct3dUSimage(sensor_data, params, c0, varargin)

% set usage defaults
num_req_input_variables = 3;
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

global Nx Ny Nt dx dy dt

% create the computational grid
Nx = params.Nx;                      % number of grid points in the x (row) direction
Ny = params.Ny;                      % number of grid points in the y (column) direction
dx = params.dx;                      % grid point spacing in the x direction [m]
dy = params.dy;                      % grid point spacing in the y direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy);

% scope parameters to make time array
Nt = params.Nt;                         % number of samples in acquisition - alternatively: size(sensor_data,3)
dt = params.dt;                         % dt between samples [s]
Nt_delay = params.trigger_delay / dt;   % number of samples to delay OR READ OUT FROM FILE NAME - HOW?
Nt = Nt_delay + Nt + params.Nt_t0_correct;
t_array = linspace(1,Nt,Nt)*dt;


%% image reconstruction

% remove the source from the time series
sensor_data = cat(3, zeros(Nx,Ny,params.Nt_zero_pad_source), sensor_data(:,:,params.Nt_zero_pad_source+1:end) );

% add zero padding for delay
sensor_data = cat(3, zeros(Nx,Ny,int32(Nt_delay)), sensor_data );

% add/remove samples from sensor_data for t0 correction
if params.Nt_t0_correct > 0
    sensor_data = cat(3, zeros(Nx,Ny,int32(params.Nt_t0_correct)), sensor_data );
elseif params.Nt_t0_correct < 0
    sensor_data = sensor_data(:,:,-params.Nt_t0_correct+1:end);
end

% zero pad the sides
if zero_pad_sides
    sensor_data = zero_padding_sides(sensor_data, Nx, Ny, dx, dy, Nt, zero_pad_sides);
end

% upsample along space (x and y)
if toUpsample
    sensor_data = upsampling_data_x2(sensor_data,Nx,Ny,dx,dy);
end

% window the data (apodising to remove artefacts due to edge of sensor)
if toApodise
    win = getWin([Nx Ny], 'Cosine');
    win = win + 0.5;
    sensor_data = bsxfun(@times, win, sensor_data);
end

% reconstruct an image using a k-space method
sensor_data = permute(sensor_data,[3 1 2]);       % reorder p_xyt to p_txy
reflection_image = kspacePlaneRecon_US(sensor_data, dx, dy, dt, c0);   % output as p_zxy
reflection_image = permute(reflection_image,[2 3 1]);               % reorder p_zxy to p_xyz

% time gain compensation
if toTimeGainCompensate
    tgc_exponent = 0.3;
    % tgc = exp(tgc_exponent * t_array * c0); % exponential
    tgc = tgc_exponent * t_array * c0; % linear
    tgc = reshape(tgc, 1, 1, length(tgc));
    reflection_image = bsxfun(@times, tgc, reflection_image);
end

% % positivity constraint (remove all negative intensities in image, e.g. back of sensor) - DOESN'T WORK SO WELL
% reflection_image(reflection_image<0) = 0;

% envelope detection per slice
if toEnvelopeDetect
    for i = 1:Ny
        reflection_image(:,i,:) = envelopeDetection(squeeze(reflection_image(:,i,:)));
    end
end

% log compression
if toLogCompress
    compression_ratio = 3;
    reflection_image = logCompression(reflection_image, compression_ratio, true);
end
 
end


function sensor_data_padded = zero_padding_sides(sensor_data, Nx, Ny, dx, dy, samples_total, pads)

    sensor_data_padded = zeros(Nx+2*pads, Ny+2*pads, samples_total);
    sensor_data_padded(pads+1:Nx+pads, pads+1:Ny+pads, :) = sensor_data;
    Nx = Nx + 2*pads;
    Ny = Ny + 2*pads;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    
end

function sensor_data_upsampled = upsampling_data_x2(sensor_data,Nx,Ny,dx,dy)

    sensor_data_upsampled = zeros( 2*Nx , 2*Ny, samples_total );
    for i = 1:Nx
        for j = 1:Ny
            sensor_data_upsampled(2*i-1,2*j-1,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i-1,2*j  ,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i  ,2*j-1,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i  ,2*j  ,:) = sensor_data(i,j,:);
        end
    end
    dx = dx/2;
    dy = dy/2;
    Nx = 2*Nx;
    Ny = 2*Ny;
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
    
