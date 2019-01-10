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

% sliceViewer


%% load set of images and find cmap suitable for each set with single bandwidth

cmins = NaN(size(bandwidths));
cmaxs = NaN(size(bandwidths));

bw_index = 1;

for bandwidth = bandwidths
    for centre_freq = centre_freqs

        reflection_image = load_image(phantom_id, centre_freq, bandwidth);
        
        meanIP = squeeze(mean(reflection_image(:,20:40,:),2));
        
        min_intensity = min(min(meanIP));
        max_intensity = max(max(meanIP));
        
        cmins(bw_index) = maybe_update_minimum( cmins(bw_index), min_intensity );
        cmaxs(bw_index) = maybe_update_maximum( cmaxs(bw_index), max_intensity );

    end
    bw_index = bw_index + 1;
end

save(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_caxes'], 'cmins', 'cmaxs')


%% replot images with same coloraxis for each set of images with single bandwidth

bw_index = 1;

for bandwidth = bandwidths
    for centre_freq = centre_freqs
        
        reflection_image = load_image(phantom_id, centre_freq, bandwidth);

        display_and_save_meanIP(reflection_image,c0,centre_freq,bandwidth,phantom_id,[cmins(bw_index),cmaxs(bw_index)])
    
    end
    bw_index = bw_index + 1;
end


%% get some SNR measures

meanIP = squeeze(mean(reflection_image(:,20:40,:),2));

meanIP_1D = mean(meanIP(40:44,1:500),1);

figure
plot(meanIP_1D)


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

function display_and_save_meanIP(reflection_image,c0,centre_freq,bandwidth,phantom_id,varargin)
% varargin to allow for optional specification of colour axis

    global kgrid t_array

    meanIP = squeeze(mean(reflection_image(:,20:40,:),2));
    
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








