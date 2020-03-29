%% data file locations & params

file_dir = '../data/imagingUS/';

% 181204 atmm with orgasol
    phantom_id = 'atmm_orgasol1_BK31[CNT]';
    file_name = '181204/atmm_orgasol1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
    trigger_delay = 0;
    samples_cut_off = 50;
    samples_t0_correct = -6;
    c0 = 1544;

% 180828 pork belly 3
%     phantom_id = 'porkBelly3_BK31[CNT]';
%     file_name = '180828\porkBelly3_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_26s20m19h_28-08-18_avg1_2D_raw.SGL';
%     trigger_delay = 0;
%     samples_cut_off = 10;
%     samples_t0_correct = -4;
%     c0 = 1460;

% 190114 lymph node (L3)
%     phantom_id = 'lymphNode2_BK31[CNT]';
%     file_name = '190114/lymphNode2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_13s51m16h_14-01-19_avg1_2D_raw.SGL';
%     trigger_delay = 0;
%     samples_cut_off = 50;
%     samples_t0_correct = -6;
%     c0 = 1520;


%% frequency bands to compound

centre_freqs = (2:1:15)*1e6;
bandwidths   = 10e6;


%% load and prepare data

[sensor_data, params] = loadSGL([file_dir file_name]);

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;

global kgrid t_array


%% generate lots of different reconstructions with varying freq filters
% and save image in .mat and meanIP in .jpg/.fig

num_images_compounded = 0;

for bandwidth = bandwidths
    for centre_freq = centre_freqs

        disp(['CENTRE FREQ ' num2str(centre_freq/1e6) ', BANDWIDTH ' num2str(bandwidth/1e6)])
        
        [reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
                                    'ZeroPad', 10, ...
                                    'Upsample', true, ...
	                                'Apodise', false, ...
                                    'FreqBandFilter', {centre_freq, bandwidth}, ...
                                    'FreqLowFilter', {}, ...
                                    'TimeGainCompensate',{}, ...
                                    'EnvelopeDetect', true, ...
                                    'LogCompress', 0, ...
                                    'SaveImageToFile', true ...
                                );
                            
        display_and_save_projection(phantom_id, reflection_image, c0, 'CentreFreq', centre_freq, 'Bandwidth', bandwidth)
        
        % add to compound or create new if first image
        if ~exist('compound_image','var')
            compound_image = reflection_image;
        else
            compound_image = compound_image + reflection_image;
        end
        num_images_compounded = num_images_compounded + 1;
        
    end
end

compound_image = compound_image / num_images_compounded;

display_and_save_projection(phantom_id,compound_image,c0)

save_compound_image(compound_image,[kgrid.dx, kgrid.dy, params.dt*c0],phantom_id)


%% compound with weights

global kgrid t_array %#ok<REDEFGG>

phantom_id = 'atmm_orgasol1_BK31[CNT]';
c0 = 1544;

centre_freqs = (2:1:15)*1e6;
bandwidths   = 10e6;

%linear weighting:
weight_2  = 0.5;
weight_15 = 1.5;
weights = linspace(weight_2,weight_15,length(centre_freqs)) / length(centre_freqs);
weighting_type = 'linear0.5-1.5';

[compound_image, voxel_size] = compound_with_weights(phantom_id,centre_freqs,bandwidths,weights);

save_compound_image(compound_image,voxel_size,[phantom_id '_weighted_' weighting_type])

% need to re-define kgrid and t_array for plotting
    Nx = size(compound_image,1);
    Ny = size(compound_image,2);
    Nt = size(compound_image,3);
    dx = voxel_size(1);
    dy = voxel_size(2);
    dt = voxel_size(3)/c0;
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    t_array = linspace(1,Nt,Nt)*dt;
display_and_save_projection(phantom_id, compound_image, c0, 'Brightness', 0.5, 'WeightingType', weighting_type)


%% LOCAL FUNCTIONS

function [reflection_image, voxel_size] = load_image(phantom_id, centre_freq, bandwidth)

    file_path = ['recon_data\' phantom_id ...
         '_f' num2str(centre_freq/1e6) ...
         '_bw' num2str(bandwidth/1e6) ...
         '.mat'];
    image_data = load(file_path);
    reflection_image = image_data.volume_data;
    voxel_size = image_data.volume_spacing;

end

function meanIP = get_mean_intensity_projection(reflection_image, xrange, yrange, zrange)

    global t_array
    
    meanIP = squeeze(mean(reflection_image(xrange,yrange,zrange),2));
    if length(t_array) > max(zrange)
        t_array = t_array(zrange);
    end

end

function maxIP = get_max_intensity_projection(reflection_image, xrange, yrange, zrange)

    global t_array
    
    maxIP = squeeze(max(reflection_image(xrange,yrange,zrange),[],2));
    if length(t_array) > max(zrange)
        t_array = t_array(zrange);
    end

end

function display_and_save_projection(phantom_id, reflection_image, c0, varargin)
% varargin to allow for optional specification of colour axis.
% compound images identified as (f,bw) = (0,0), saved as '_compound'
% for any other images (f,bw) is specified in varargin

    global kgrid t_array
    
    num_req_input_variables = 3;
    centre_freq = 0;
    bandwidth = 0;
    
    % check for optional inputs
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'CentreFreq'
                    centre_freq = varargin{input_index+1};
                case 'Bandwidth'
                    bandwidth = varargin{input_index+1};
                case 'Colormap'
                    climits = varargin{input_index+1};
                case 'Brightness'
                    brightness = varargin{input_index+1};
                case 'WeightingType'
                    weighting_type = varargin{input_index+1};
                otherwise
                    error('Unknown optional input.');
            end
        end
    end
    
    % define limits for intensity projections
    switch phantom_id
        case 'atmm_orgasol1_BK31[CNT]'
            xrange = 10:210;
            yrange = 70:110;
            zrange = 105:753;
            mIP = get_mean_intensity_projection(reflection_image,xrange,yrange,zrange);
        case 'porkBelly3_BK31[CNT]'
            xrange = 22:322;
            yrange = 145:185;
            zrange = 1:377;
            mIP = get_mean_intensity_projection(reflection_image,xrange,yrange,zrange);
        case 'lymphNode2_BK31[CNT]'
            xrange = 22:322;
            yrange = 160:200; % 5th slice % full volume = 21:321
            zrange = 30:360;
            mIP = get_max_intensity_projection(reflection_image,xrange,yrange,zrange);
    end
	
    % plot figure
    figure
    imagesc(kgrid.x_vec(xrange)*1e3,t_array*c0*1e3,mIP')
        axis image
        colormap(gray)
        if exist('climits','var')
            caxis(climits)
        end
        if exist('brightness','var')
            brighten(brightness)
        end
        % title
        if centre_freq == 0 && bandwidth == 0
            title('z-x MeanIP 2 mm slice compound')
            if exist('weighting_type','var')
                title(['z-x MeanIP 2 mm slice compound weighted ' weighting_type])
            end
        else
            switch phantom_id
                case 'lymphNode2_BK31[CNT]'
                    title(['z-x MaxIP 2 mm slice: f = ' num2str(centre_freq/1e6) ', bw = ' num2str(bandwidth/1e6)])
                otherwise
                    title(['z-x MeanIP 2 mm slice: f = ' num2str(centre_freq/1e6) ', bw = ' num2str(bandwidth/1e6)])
            end
        end
        xlabel('x-position [mm]')
        ylabel('z-position [mm]')
        set(gca,'FontName','Arial')
        set(gca,'FontSize',12)
    
	% save figure
    file_img = ['..\figures\_Matlab figs\freqCompounding\' phantom_id];
    if centre_freq == 0 && bandwidth == 0
        if exist('weighting_type','var')
            file_img = [file_img '_weighted_' weighting_type];
        end
        file_img = [file_img '_compound'];
    else
        file_img = [file_img ...
                     '_f' num2str(centre_freq/1e6) ...
                     '_bw' num2str(bandwidth/1e6)   ];
    end
    saveas(gcf,[file_img '.fig'])
    saveas(gcf,[file_img '.jpg'])

end

function save_compound_image(volume_data, volume_spacing, phantom_id)

    tic
    disp('Saving compound image data ..')
    
    file_path = ['recon_data\' phantom_id '_compound.mat'];
    save(file_path,'volume_data','volume_spacing','-v7.3')
    
    disp(['  completed in ' scaleTime(toc)]);

end

function [compound_image, voxel_size] = compound_with_weights(phantom_id, centre_freqs, bandwidths, weights)
% currently only supports weights for different centre_freqs not for
% different bandwidths

    disp(['Compounding ' phantom_id ' with weights:'])
    
    assert(sum(weights)==1);
    
    for bandwidth = bandwidths
        for idx_f = 1:length(centre_freqs)

            centre_freq = centre_freqs(idx_f);
            weight = weights(idx_f);

            tic
            disp(['Loading and adding f = ' num2str(centre_freq/1e6) ', bw = ' num2str(bandwidth/1e6) ' MHz ..'])

            % load frequency banded image
            [reflection_image, voxel_size] = load_image(phantom_id, centre_freq, bandwidth);
            % add image with weight
            if ~exist('compound_image','var')
                compound_image = weight * reflection_image;
            else
                compound_image = compound_image + weight * reflection_image;
            end
            
            disp(['  completed in ' scaleTime(toc)])
            
        end
    end
    
end

