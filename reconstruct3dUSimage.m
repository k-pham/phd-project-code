function [reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, varargin)

%% read from varargin

% set usage defaults
num_req_input_variables = 3;
zero_pad_sides = 0;
toUpsample = false;
toApodise = false;
freq_band = [];
tgc_params = {};
toEnvelopeDetect = false;
toLogCompress = 0;

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
            case 'FreqBandFilter'
                freq_band = varargin{input_index + 1};
                assert(length(freq_band)==2,'Need list with 2 cut-off frequencies.')
            case 'TimeGainCompensate'
                tgc_params = varargin{input_index + 1};
                assert(iscell(tgc_params),'Need cell with {method, strength}.')
            case 'EnvelopeDetect'
                toEnvelopeDetect = varargin{input_index + 1};
            case 'LogCompress'
                toLogCompress = varargin{input_index + 1};
                assert(isnumeric(toLogCompress),'Need number for compression ratio.')
            otherwise
                error('Unknown optional input.');
        end
    end
end


%% assign parameters from params to use in script and globally

global Nx Ny Nt dx dy dt kgrid t_array

% create the computational grid
Nx = params.Nx;                      % number of grid points in the x (row) direction
Ny = params.Ny;                      % number of grid points in the y (column) direction
dx = params.dx;                      % grid point spacing in the x direction [m]
dy = params.dy;                      % grid point spacing in the y direction [m]

% scope parameters to make time array
Nt = size(sensor_data,3);               % number of samples in acquisition - alternatively: params.Nt
dt = params.dt;                         % dt between samples [s]
Nt_delay = params.trigger_delay / dt;   % number of samples to delay OR READ OUT FROM FILE NAME - HOW?


%% prepare data for reconstruction

% zero pad source noise in the time series
if params.Nt_zero_pad_source
    sensor_data = zero_padding_source(sensor_data, params.Nt_zero_pad_source);
end

% add zero padding for delay
if Nt_delay ~= 0
    sensor_data = zero_padding_delay(sensor_data, Nt_delay);
end

% add/remove samples from sensor_data for t0 correction
if params.Nt_t0_correct
    sensor_data = correcting_t0(sensor_data, params.Nt_t0_correct);
end

% zero pad the sides
if zero_pad_sides
    sensor_data = zero_padding_sides(sensor_data, zero_pad_sides);
end

% upsample along space (x and y)
if toUpsample
    sensor_data = upsampling_x2(sensor_data);
end

% apodising data (to remove edge wave artefacts)
if toApodise
    sensor_data = apodising(sensor_data);
end

% frequency band filtering data (to use for frequency compounding)
if ~isempty(freq_band)
    sensor_data = freq_filtering(sensor_data, freq_band);
end


%% reconstruct image using k-space method

sensor_data = permute(sensor_data,[3 1 2]);                             % reorder p_xyt to p_txy
reflection_image = kspacePlaneRecon_US(sensor_data, dx, dy, dt, c0);    % output as p_zxy
reflection_image = permute(reflection_image,[2 3 1]);                   % reorder p_zxy to p_xyz


%% update kgrid and make t_array for use outside

kgrid = kWaveGrid(Nx, dx, Ny, dy);
t_array = linspace(1,Nt,Nt)*dt;


%% image processing steps

if ~isempty(tgc_params)
    reflection_image = time_gain_compensating(reflection_image, c0, tgc_params);
end

if toEnvelopeDetect
    reflection_image = envelope_detecting(reflection_image);
end

if toLogCompress
    reflection_image = log_compressing(reflection_image, toLogCompress);
end

end


%% zero padding source
function sensor_data_padded = zero_padding_source(sensor_data, pads)

    global Nx Ny
    
    sensor_data_padded = cat(3, zeros(Nx,Ny,pads), sensor_data(:,:,pads+1:end) );
    assert(isequal( size(sensor_data_padded), size(sensor_data) ))

end


%% zero padding delay
function sensor_data_padded = zero_padding_delay(sensor_data, pads)

    global Nx Ny Nt
    
    sensor_data_padded = cat(3, zeros(Nx,Ny,int32(pads)), sensor_data );
    assert(isequal( size(sensor_data_padded), size(sensor_data)+[0,0,pads] ))
    
    Nt = Nt + pads;

end


%% correcting t0
function sensor_data_t0_corrected = correcting_t0(sensor_data, correction)

    global Nx Ny Nt
    
    if correction > 0
        sensor_data_t0_corrected = cat(3, zeros(Nx,Ny,correction), sensor_data);
    elseif correction < 0
        sensor_data_t0_corrected = sensor_data(:,:,-correction+1:end);
    end
	assert(isequal( size(sensor_data_t0_corrected), size(sensor_data)+[0,0,correction] ))
    
    Nt = Nt + correction;

end


%% zero padding sides
function sensor_data_padded = zero_padding_sides(sensor_data, pads)

    global Nx Ny Nt
    
    sensor_data_padded = zeros(Nx+2*pads, Ny+2*pads, Nt);
    sensor_data_padded(pads+1:Nx+pads, pads+1:Ny+pads, :) = sensor_data;
    assert(isequal( size(sensor_data_padded), size(sensor_data)+[2*pads,2*pads,0] ))
    
    Nx = Nx + 2*pads;
    Ny = Ny + 2*pads;
    
end


%% upsampling data *2
function sensor_data_upsampled = upsampling_x2(sensor_data)

    global Nx Ny Nt dx dy

    sensor_data_upsampled = zeros( 2*Nx, 2*Ny, Nt );
    for i = 1:Nx
        for j = 1:Ny
            sensor_data_upsampled(2*i-1,2*j-1,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i-1,2*j  ,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i  ,2*j-1,:) = sensor_data(i,j,:);
            sensor_data_upsampled(2*i  ,2*j  ,:) = sensor_data(i,j,:);
        end
    end
    assert(isequal( size(sensor_data_upsampled), [2*Nx, 2*Ny, Nt] ))
    
	Nx = 2*Nx;
    Ny = 2*Ny;
    dx = dx/2;
    dy = dy/2;

end


%% apodising data
function sensor_data_apodised = apodising(sensor_data)

    global Nx Ny
    
    win = getWin([Nx Ny], 'Cosine');
    win = win + 0.5;
    sensor_data_apodised = bsxfun(@times, win, sensor_data);
    assert(isequal( size(sensor_data_apodised), size(sensor_data) ))
    
end


%% time gain compensation
function reflection_image_tgc = time_gain_compensating(reflection_image, c0, tgc_params)

    global t_array
    
    disp('Time gain compensating ...')
    tic
    
    [method, strength] = tgc_params{:};
    switch method
        case 'Linear'
            tgc = strength * t_array * c0;
        case 'Exponential'
            tgc = exp(strength * t_array * c0);
    end    
    tgc = reshape(tgc, 1, 1, length(tgc));
    reflection_image_tgc = bsxfun(@times, tgc, reflection_image);
    assert(isequal( size(reflection_image_tgc), size(reflection_image) ))
    
    disp(['  completed in ' scaleTime(toc)]);
    
end


%% envelope detection
function reflection_image_env = envelope_detecting(reflection_image)

    global Ny
    
    disp('Envelope detecting ...')
    tic
    
    reflection_image_env = zeros(size(reflection_image));
    for i = 1:Ny
        reflection_image_env(:,i,:) = envelopeDetection(squeeze(reflection_image(:,i,:)));
    end
    assert(isequal( size(reflection_image_env), size(reflection_image) ))
    
    disp(['  completed in ' scaleTime(toc)]);

end


%% log compression
function reflection_image_log = log_compressing(reflection_image, compression_ratio)

    disp('Log Compressing ...')
    tic
    
    reflection_image_log = logCompression(reflection_image, compression_ratio, true);
    assert(isequal( size(reflection_image_log), size(reflection_image) ))
    
    disp(['  completed in ' scaleTime(toc)]);

end


%% frequency bandpass filtering data (for frequency compounding)
function sensor_data_filtered = freq_filtering(sensor_data, freq_band)

    global Nx Ny dt
    
    disp('Frequency bandpass filtering ...'),
    tic
    
    for i = 1:Nx
        for j = 1:Ny
            sensor_data_filtered(i,j,:) = applyFilter(squeeze(sensor_data(i,j,:)),1/dt,freq_band,'BandPass');
        end
    end
    
    disp(['  completed in ' scaleTime(toc)]);

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
    
