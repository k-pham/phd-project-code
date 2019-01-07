function [reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, varargin)

%% read from varargin

% set usage defaults
num_req_input_variables = 3;
zero_pad_sides = 0;
toUpsample = false;
toApodise = false;
freqfilter_params = {};
tgc_params = {};
toEnvelopeDetect = false;
compression_ratio = 0;
toSaveImage = false;

% replace with user defined values if provided
if nargin < num_req_input_variables
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'ZeroPad'
                zero_pad_sides = varargin{input_index + 1};
                assert(isnumeric(zero_pad_sides),'Need number for number of pads.')
            case 'Upsample'
                toUpsample = varargin{input_index + 1};
            case 'Apodise'
                toApodise = varargin{input_index + 1};
            case 'FreqBandFilter'
                freqfilter_params = varargin{input_index + 1};
                assert(iscell(freqfilter_params),'Need cell array {centre_freq(Hz) bandwidth(Hz)}.')
                assert(length(freqfilter_params)==2,'Need cell array of length 2.')
            case 'TimeGainCompensate'
                tgc_params = varargin{input_index + 1};
                assert(iscell(tgc_params),'Need cell array {method, strength}.')
                assert(length(tgc_params)==2,'Need cell array of length 2.')
            case 'EnvelopeDetect'
                toEnvelopeDetect = varargin{input_index + 1};
            case 'LogCompress'
                compression_ratio = varargin{input_index + 1};
                assert(isnumeric(compression_ratio),'Need number for compression ratio.')
            case 'SaveImageToFile'
                toSaveImage = varargin{input_index + 1};
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
if ~isempty(freqfilter_params)
    sensor_data = freq_filtering(sensor_data, freqfilter_params);
end


%% reconstruct image using k-space method

sensor_data = permute(sensor_data,[3 1 2]);                             % reorder p_xyt to p_txy
reflection_image = kspacePlaneRecon_US(sensor_data, dx, dy, dt, c0);    % output as p_zxy
reflection_image = permute(reflection_image,[2 3 1]);                   % reorder p_zxy to p_xyz


%% image trimming & interpolation

reflection_image = trim_image_x(reflection_image,xStart,xEnd);
reflection_image = trim_image_y(reflection_image,yStart,yEnd);
reflection_image = trim_image_z(reflection_image,zStart,zEnd);


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

if compression_ratio ~= 0
    reflection_image = log_compressing(reflection_image, compression_ratio);
end


%% saving image data to .mat

if toSaveImage
    savingImageToMat(reflection_image, params.file_data, dt*c0, freqfilter_params)     % omit factor 1/2 in dz because of doubled depth bug
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
function sensor_data_filtered = freq_filtering(sensor_data, freqfilter_params)

    global Nx dt
    
    disp('Frequency bandpass filtering ...'),
    tic
    
    [centre_freq, bandwidth] = freqfilter_params{:};
    bandwidth_pc = bandwidth / centre_freq * 100;       % convert bandwidth to % of centre frequency for gaussianFilter
    sensor_data_filtered = zeros(size(sensor_data));
    
    for i = 1:Nx
        sensor_data_filtered(i,:,:) = gaussianFilter(squeeze(sensor_data(i,:,:)),1/dt,centre_freq,bandwidth_pc);
    end
    
    disp(['  completed in ' scaleTime(toc)]);

end


%% saving image data to .mat (for sliceViewer)
function savingImageToMat(reflection_image, file_data, dz, freqfilter_params)

    global Nx Ny kgrid
    
    disp('Saving data to .mat ...'),
    tic

    Nz = size(reflection_image,3);
    volume_data = reshape(reflection_image,Nx,Ny,Nz);       %#ok<NASGU>
    volume_spacing = [kgrid.dx, kgrid.dy, dz];              %#ok<NASGU>

    phantom_id = strtok(file_data,'@'); % parse string up to specified delimiter
    phantom_id = phantom_id(8:end);     % remove date folder from string
    [centre_freq, bandwidth] = freqfilter_params{:};
    file_image = ['recon_data\' phantom_id ...
                    '_f' num2str(centre_freq/1e6) ... % incl freq filtering info ...
                    '_bw' num2str(bandwidth/1e6) ...  % in image file name
                    '.mat'];
    
    save(file_image,'volume_data','volume_spacing','-v7.3')
    
    disp(['  completed in ' scaleTime(toc)]);
    
end


%% trim image in x direction
function reflection_image_trimmed = trim_image_x(reflection_image,xStart,xEnd)

    global Nx

    reflection_image_trimmed = reflection_image(xStart:xEnd,:,:);
    
    Nx = size(reflection_image_trimmed,1);

end


%% trim image in y direction
function reflection_image_trimmed = trim_image_y(reflection_image,yStart,yEnd)

    global Ny

    reflection_image_trimmed = reflection_image(:,yStart:yEnd,:);
    
    Ny = size(reflection_image_trimmed,2);

end


%% trim image in z direction
function reflection_image_trimmed = trim_image_z(reflection_image,zStart,zEnd)

    global Nz

    reflection_image_trimmed = reflection_image(:,:,zStart:zEnd);
    
    Nz = size(reflection_image_trimmed,3);

end











