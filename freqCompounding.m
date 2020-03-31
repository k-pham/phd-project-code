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
bandwidths   = 2e6;


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
                            
        display_and_save_projection(phantom_id, reflection_image, c0, 'CentreFreq', centre_freq, 'Bandwidth', bandwidth, 'Brightness', 0.5)
        
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

display_and_save_projection(phantom_id, compound_image, c0, 'Brightness', 0.5)

save_compound_image(phantom_id, compound_image, [kgrid.dx, kgrid.dy, params.dt*c0])


%% compound with weights

global kgrid t_array %#ok<REDEFGG>

phantom_id = 'atmm_orgasol1_BK31[CNT]';
c0 = 1544;

centre_freqs = (2:1:15)*1e6;
bandwidths   = 2e6;

% linear weighting:
weight_2  = 0.5;
weight_15 = 1.5;
weights = linspace(weight_2,weight_15,length(centre_freqs)) / length(centre_freqs);
weighting_type = ['bw' num2str(bandwidths/1e6) '_linear' num2str(weight_2) '-' num2str(weight_15)];

% compound with weights and save image data
[compound_image, voxel_size] = compound_with_weights(phantom_id,centre_freqs,bandwidths,weights);
save_compound_image([phantom_id '_weighted_' weighting_type], compound_image, voxel_size)

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

% plot showing weights and save
figure
bar(centre_freqs/1e6,weights)
    xlim([1,16])
    ylim([0,0.12])
    xlabel('centre frequency / MHz')
    ylabel('weight')
    set(gca,'FontSize',14)
    save(gcf,['..\figures\_Matlab figs\freqCompounding\weight_linear' num2str(weight_2) '-' num2str(weight_15) '.jpg'])


%% assess image quality of weighted compounds

% file id to get weighting_type
phantom_id = 'atmm_orgasol1_BK31[CNT]';
bandwidths   = 10e6;
weight_2  = 1;
weight_15 = 1;
weighting_type = ['bw' num2str(bandwidths/1e6) '_linear' num2str(weight_2) '-' num2str(weight_15) '_compound.mat'];

% load image data with specified weighting_type
file_path = ['recon_data\' phantom_id '_weighted_' weighting_type];
image_data = load(file_path);
compound_image = image_data.volume_data;
voxel_size = image_data.volume_spacing;

% get meanIP
xrange = 10:210;
yrange = 70:110;
zrange = 105:753;
meanIP = get_mean_intensity_projection(compound_image, xrange, yrange, zrange);
figure
imagesc(meanIP')

% extract resolution/signal/scattering from meanIP
[resoLat, resoAxi]                      = get_resolution_of_tube(meanIP, voxel_size(1), voxel_size(3));
signal_tube                             = get_peak_signal_of_tube(meanIP);
[scatter_water_mean, scatter_water_std] = get_scattering_distr_in_water(meanIP);
[scatter_atmm_mean, scatter_atmm_std]   = get_scattering_distr_in_atmm(meanIP);

% calculate image quality metrics
specSNR = signal_tube       / scatter_water_mean;
scatSNR = scatter_atmm_mean / scatter_water_mean;
scatCNR = (scatter_atmm_mean - scatter_water_mean) / (scatter_atmm_std + scatter_water_std);


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

function save_compound_image(phantom_id, volume_data, volume_spacing)

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

function [resoLat, resoAxi] = get_resolution_of_tube(mIP, dx, dz)

    ROI = mIP(106:130,300:370);
    
    [~,peak_pos] = max(ROI,[],'all','linear');
    [peak_pos_x, peak_pos_z] = ind2sub(size(ROI),peak_pos);
    
    profile_x = ROI(:,peak_pos_z);
    profile_z = ROI(peak_pos_x,:);
    
    resoLat = fwhm(profile_x,dx,1);
    resoAxi = fwhm(profile_z,dz,1);

end

function signal_tube = get_peak_signal_of_tube(mIP) % hard numbers !
% limits customised to third tube

    signal_tube_front = max(max(mIP(113:123,300:370)));
    signal_tube_back  = max(max(mIP(113:123,370:440)));
    signal_tube = mean([signal_tube_front,signal_tube_back]);

end

function [scatter_water_mean, scatter_water_std] = get_scattering_distr_in_water(mIP) % hard numbers !
% limits customised to third tube

    ROI = mIP(113:123,355:395);
    
    scatter_water_mean = mean(ROI(:));
    scatter_water_std  = std(ROI(:));

end

function [scatter_atmm_mean, scatter_atmm_std] = get_scattering_distr_in_atmm(mIP) % hard numbers !
% limits customised to third tube

    ROI_left  = mIP(90:100,355:395);
    ROI_right = mIP(136:146,355:395);
    ROI_front = mIP(113:123,270:310);
    ROI_back  = mIP(113:123,440:480);
    
    ROI = [ROI_left ROI_right ROI_front ROI_back];

    scatter_atmm_mean = mean(ROI(:));
    scatter_atmm_std  = std(ROI(:));

end
