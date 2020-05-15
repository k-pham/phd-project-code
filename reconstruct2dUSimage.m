function [reflection_image] = reconstruct2dUSimage(sensor_data, params, c0, varargin)

%% read from varargin

% set usage defaults
num_req_input_variables = 3;
toUpsample = false;
toApodise = false;
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
            case 'Upsample'
                toUpsample = varargin{input_index + 1};
            case 'Apodise'
                toApodise = varargin{input_index + 1};
            case 'TimeGainCompensate'
                tgc_params = varargin{input_index + 1};
                assert(iscell(tgc_params),'Need cell array {method, strength}.')
                if ~isempty(tgc_params)
                    assert(length(tgc_params)==2,'Need cell array of length 2.')
                end
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


%% prepare data for reconstruction

% zero pad source noise in the time series
if params.Nt_zero_pad_source
    sensor_data = [zeros(Nx,params.Nt_zero_pad_source) sensor_data(:,params.Nt_zero_pad_source+1:end)];
end

% add zero padding for delay
if Nt_delay ~= 0
    sensor_data = [zeros(Nx,int32(Nt_delay)) sensor_data];
    Nt = Nt + Nt_delay;
end

% add/remove samples from sensor_data for t0 correction
if params.Nt_t0_correct > 0
    sensor_data = [ zeros(Nx,int32(params.Nt_t0_correct)) sensor_data ];
elseif params.Nt_t0_correct < 0
    sensor_data = sensor_data(:,-params.Nt_t0_correct+1:end);
end
Nt = Nt + params.Nt_t0_correct;

% upsample along space (x)
if toUpsample
    sensor_data_upsampled = zeros( 2*Nx , Nt );
    for i = 1:Nx
        sensor_data_upsampled(2*i-1,:) = sensor_data(i,:);
        sensor_data_upsampled(2*i  ,:) = sensor_data(i,:);
    end
    sensor_data = sensor_data_upsampled;
    clear sensor_data_upsampled
    dx = dx/2;
    Nx = 2*Nx;
end

% apodising data (to remove edge wave artefacts)
if toApodise
    win = getWin(Nx, 'Cosine');
    sensor_data = bsxfun(@times, win, sensor_data);
end


%% reconstruct an image using a k-space method

reflection_image = kspaceLineRecon_US(sensor_data', dx, dt, c0);            % input p_tx, output p_zx
reflection_image = permute(reflection_image,[2 1]);                         % reorder p_zx to p_xz


%% trim image in z by half

reflection_image = reflection_image(:,1:round(Nt/2));
Nt = size(reflection_image,2);


%% update kgrid and make t_array for use outside

kgrid = kWaveGrid(Nx, dx, Ny, dy);
t_array = linspace(1,Nt,Nt)*dt;


%% image processing steps

% time gain compensation
if ~isempty(tgc_params)
    [method, strength] = tgc_params{:};
    switch method
        case 'Linear'
            tgc = strength * t_array * c0;
        case 'Exponential'
            tgc = exp(strength * t_array * c0);
    end
    tgc = reshape(tgc, 1, length(tgc));
    reflection_image = bsxfun(@times, tgc, reflection_image);
end

% envelope detection
if toEnvelopeDetect
    reflection_image = envelopeDetection(squeeze(reflection_image));            % p_xz
end

% log compression
if compression_ratio ~= 0
    reflection_image = logCompression(reflection_image, compression_ratio, true);   % p_xz
end


%% saving image data to .mat

if toSaveImage
    [Nz] = size(reflection_image,2);
    volume_data = reshape(reflection_image,Nx,1,Nz);
    volume_spacing = [kgrid.dx, kgrid.dx, dt*c0];           % omit factor 1/2 in dz because of doubled depth bug

    phantom_id = strtok(params.file_data,'@');  % parse string up to specified delimiter
    phantom_id = phantom_id(8:end);             % remove date folder from string
    save(['recon_data\' phantom_id '.mat'],'volume_data','volume_spacing','-v7.3')
end


end
