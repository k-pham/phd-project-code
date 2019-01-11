%% loading sensor data & preparing for reconstruction

clear all %#ok<CLALL>
run('USimagingPhantoms.m')

phantom_id = 'atmm_orgasol1_BK31[CNT]';

[sensor_data, params] = loadSGL([file_dir file_name]);

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;

global Nx Ny kgrid t_array %#ok<NUSED>

bandwidths   = [1:1:3,4:2:10,15:5:40] *1e6;
centre_freqs = (1:1:35) *1e6;


%% generate lots of different reconstructions with varying freq filters
% and save image in .mat and meanIP in .jpg/.fig

tic

for bandwidth = bandwidths
    for centre_freq = centre_freqs

        [reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
                                    'ZeroPad', 10, ...
                                    'FreqBandFilter', {centre_freq, bandwidth}, ...
                                    'TimeGainCompensate',{'Exponential',50}, ...
                                    'EnvelopeDetect', true, ...
                                    'SaveImageToFile', true ...
                                );

        display_and_save_meanIP(reflection_image, c0, centre_freq, bandwidth, phantom_id)

    end
end

disp(['  completed in ' scaleTime(toc)]);

% sliceViewer


%% load set of images and find cmap suitable for each set with single bandwidth

tic

cmins = NaN(size(bandwidths));
cmaxs = NaN(size(bandwidths));

bw_index = 1;
for bandwidth = bandwidths
    for centre_freq = centre_freqs

        reflection_image = load_image(phantom_id, centre_freq, bandwidth);
        
        meanIP = get_mean_intensity_projection(reflection_image);
        
        min_intensity = min(min(meanIP));
        max_intensity = max(max(meanIP));
        
        cmins(bw_index) = maybe_update_minimum( cmins(bw_index), min_intensity );
        cmaxs(bw_index) = maybe_update_maximum( cmaxs(bw_index), max_intensity );

    end
    bw_index = bw_index + 1;
end

save(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_caxes'], 'cmins', 'cmaxs')

disp(['  completed in ' scaleTime(toc)]);


%% replot images with same coloraxis for each set of images with single bandwidth

tic

bw_index = 1;
for bandwidth = bandwidths
    for centre_freq = centre_freqs
        
        reflection_image = load_image(phantom_id, centre_freq, bandwidth);

        display_and_save_meanIP(reflection_image,c0,centre_freq,bandwidth,phantom_id,[cmins(bw_index),cmaxs(bw_index)])
    
    end
    bw_index = bw_index + 1;
end

disp(['  completed in ' scaleTime(toc)]);


%% estimate SNR and CNR measures

tic

signal_tube_ar   = zeros(length(bandwidths),length(centre_freqs));
scatter_water_ar = zeros(length(bandwidths),length(centre_freqs));
scatter_atmm_ar  = zeros(length(bandwidths),length(centre_freqs));

bw_index = 1;
for bandwidth = bandwidths

    cf_index = 1;
    for centre_freq = centre_freqs

        reflection_image = load_image(phantom_id, centre_freq, bandwidth);

        meanIP = get_mean_intensity_projection(reflection_image);
        %figure
        %imagesc(meanIP')

        %show_axial_and_lateral_profiles(meanIP)

        signal_tube   = get_peak_signal_of_tube(meanIP);
        scatter_water = get_avg_scattering_in_water(meanIP);
        scatter_atmm  = get_avg_scattering_in_atmm(meanIP);

        signal_tube_ar(bw_index,cf_index)   = signal_tube;
        scatter_water_ar(bw_index,cf_index) = scatter_water;
        scatter_atmm_ar(bw_index,cf_index)  = scatter_atmm;
        
        %SNR = signal_tube / scatter_water;
        %CNR = scatter_atmm / scatter_water;

        cf_index = cf_index + 1;
    end
    bw_index = bw_index + 1;
end

save(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_imageAnalysis'], 'signal_tube_ar', 'scatter_water_ar', 'scatter_atmm_ar')

SNR = signal_tube_ar ./ scatter_water_ar;
CNR = scatter_atmm_ar ./ scatter_water_ar;

save(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_imageQuality'], 'SNR', 'CNR')

disp(['  completed in ' scaleTime(toc)]);


%% compounding

compound_image = zeros([Nx,Ny,Nz]);

for bandwidth = bandwidths
    for centre_freq = centre_freqs

        reflection_image = load_image(phantom_id, centre_freq, bandwidth);

        compound_image = compound_image + reflection_image;

    end
    compound_image = compound_image / length(centre_freq);
end

volume_data = reshape(compound_image,Nx,Ny,Nz);
volume_spacing = [kgrid.dx, kgrid.dy, params.dt*c0];

file_path = ['recon_data\' phantom_id '_compound.mat'];
save(file_path,'volume_data','volume_spacing','-v7.3')

sliceViewer


%% local functions

function reflection_image = load_image(phantom_id, centre_freq, bandwidth)

    file_path = ['recon_data\frequency compounding\' phantom_id ...
         '_f' num2str(centre_freq/1e6) ...
         '_bw' num2str(bandwidth/1e6) ...
         '.mat'];
    image_data = load(file_path);
    reflection_image = image_data.volume_data;

end

function meanIP = get_mean_intensity_projection(reflection_image) % hard numbers !
% limits customised to atmm_orgasol1_BK31[CNT]

    meanIP = squeeze(mean(reflection_image(:,20:40,:),2));

end

function display_and_save_meanIP(reflection_image,c0,centre_freq,bandwidth,phantom_id,varargin)
% varargin to allow for optional specification of colour axis

    global kgrid t_array

    meanIP = get_mean_intensity_projection(reflection_image);
    
    figure
    imagesc(kgrid.x_vec*1e3,t_array*c0*1e3,meanIP')
        axis image
        colormap(gray)
        if ~isempty(varargin{1})
            caxis(varargin{1})
        end
        brighten(0.5)
        title(['z-x MeanIP 1.5 mm slice: f = ' num2str(centre_freq/1e6) ', bw = ' num2str(bandwidth/1e6)])
        xlabel('x-position [mm]')
        ylabel('z-position [mm]')
        set(gca,'FontName','Arial')
        set(gca,'FontSize',12)

    file_img = ['..\figures\_Matlab figs\freqCompounding\' phantom_id ...
         '_f' num2str(centre_freq/1e6) ...
         '_bw' num2str(bandwidth/1e6)   ];
    saveas(gcf,[file_img '.fig'])
    saveas(gcf,[file_img '.jpg'])

end

function min = maybe_update_minimum(min, new_value)

    if new_value < min || isnan(min)
        min = new_value;
    end

end

function max = maybe_update_maximum(max, new_value)

    if new_value > max || isnan(max)
        max = new_value;
    end

end

function show_axial_and_lateral_profiles(meanIP) % hard numbers !
    
    % second tube
    % meanIP_1D_axi = mean(meanIP(40:44,:),1);
    % meanIP_1D_lat = mean(meanIP(:,225:265),2);

    % third tube
    meanIP_1D_axi = mean(meanIP(61:65,:),1);
    meanIP_1D_lat = mean(meanIP(:,385:425),2);

    figure
    plot(meanIP_1D_axi)
    figure
    plot(meanIP_1D_lat)
    
end

function signal_tube = get_peak_signal_of_tube(meanIP) % hard numbers !
% limits customised to third tube

    signal_tube_front = max(max(meanIP(61:65,300:400)));
    signal_tube_back  = max(max(meanIP(61:65,400:500)));
    signal_tube = mean([signal_tube_front,signal_tube_back]);

end

function scatter_water = get_avg_scattering_in_water(meanIP) % hard numbers !
% limits customised to third tube

    scatter_water = mean(mean(meanIP(61:65,385:425)));

end

function scatter_atmm = get_avg_scattering_in_atmm(meanIP) % hard numbers !
% limits customised to third tube

    scatter_atmm_left  = mean(mean(meanIP(51:55,385:425)));
    scatter_atmm_right = mean(mean(meanIP(72:76,385:425)));
    scatter_atmm_front = mean(mean(meanIP(61:65,300:340)));
    scatter_atmm_back  = mean(mean(meanIP(61:65,465:505)));
    scatter_atmm = mean([scatter_atmm_left,scatter_atmm_right,scatter_atmm_front,scatter_atmm_back]);

end


