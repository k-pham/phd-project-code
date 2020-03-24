%% data file locations & params

% 181204 agar-based tissue mimicking phantom (atmm) with orgasol - BK31[CNT]
phantom_id = 'atmm_orgasol1_BK31[CNT]';
file_name = '181204/atmm_orgasol1_BK31[CNT]@0nm_t[0]_dx[100µm]_dy[100µm]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
c0 = 1544;
trigger_delay = 0;
samples_cut_off = 50;
samples_t0_correct = -6;


%% frequency bands to compound

bandwidths   = [1:1:3,4:2:10,15:5:40] *1e6;
centre_freqs = (1:1:35) *1e6;


%% load and prepare data

[sensor_data, params] = loadSGL([file_dir file_name]);

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;

global Nx Ny Nt kgrid t_array %#ok<NUSED>


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

sliceViewer


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

signal_tube_ar        = zeros(length(bandwidths),length(centre_freqs));
scatter_water_mean_ar = zeros(length(bandwidths),length(centre_freqs));
scatter_water_std_ar  = zeros(length(bandwidths),length(centre_freqs));
scatter_atmm_mean_ar  = zeros(length(bandwidths),length(centre_freqs));
scatter_atmm_std_ar   = zeros(length(bandwidths),length(centre_freqs));

bw_index = 1;
for bandwidth = bandwidths

    cf_index = 1;
    for centre_freq = centre_freqs

        reflection_image = load_image(phantom_id, centre_freq, bandwidth);

        meanIP = get_mean_intensity_projection(reflection_image);
        figure
        imagesc(meanIP')

        show_axial_and_lateral_profiles(meanIP)

        signal_tube                             = get_peak_signal_of_tube(meanIP);
        [scatter_water_mean, scatter_water_std] = get_scattering_distr_in_water(meanIP);
        [scatter_atmm_mean, scatter_atmm_std]   = get_scattering_distr_in_atmm(meanIP);

        signal_tube_ar(bw_index,cf_index)        = signal_tube;
        scatter_water_mean_ar(bw_index,cf_index) = scatter_water_mean;
        scatter_water_std_ar(bw_index,cf_index)  = scatter_water_std;
        scatter_atmm_mean_ar(bw_index,cf_index)  = scatter_atmm_mean;
        scatter_atmm_std_ar(bw_index,cf_index)   = scatter_atmm_std;
        
        cf_index = cf_index + 1;
    end
    bw_index = bw_index + 1;
end

save(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_imageAnalysis'], ...
    'signal_tube_ar', 'scatter_water_mean_ar', 'scatter_water_std_ar', 'scatter_atmm_mean_ar', 'scatter_atmm_std_ar')

specSNR = signal_tube_ar       ./ scatter_water_mean_ar;
scatSNR = scatter_atmm_mean_ar ./ scatter_water_mean_ar;
CNR = (scatter_atmm_mean_ar - scatter_water_mean_ar) ./ (scatter_atmm_std_ar + scatter_water_std_ar);

save(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_imageQuality'], 'specSNR', 'scatSNR', 'CNR')

disp(['  completed in ' scaleTime(toc)]);


%% plot SNR and CNR measures in 2D vs centre freq and bandwidth

load(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_imageAnalysis'])
load(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_imageQuality'])

plot_fields = {signal_tube_ar, scatter_water_mean_ar, scatter_water_std_ar, scatter_atmm_mean_ar, scatter_atmm_std_ar, ...
                specSNR, ...
                scatSNR, ...
                CNR};
plot_titles = {'signal tube', 'scatter water mean', 'scatter water std', 'scatter atmm mean', 'scatter atmm std', ...
                'specSNR = signal tube / scatter water mean', ...
                'scatSNR = scatter atmm mean / scatter water mean', ...
                'CNR = (atmm mean - water mean) / (atmm std + water std)'};
img_files = {'signal_tube', 'scatter_water_mean', 'scatter_water_std', 'scatter_atmm mean', 'scatter_atmm_std', ...
                'specSNR', ...
                'scatSNR', ...
                'CNR'};

display_and_save_multiple_fields_vs_bw_cf(plot_fields, plot_titles, img_files)


%% select which images to use for compounding

load(['..\figures\_Matlab figs\freqCompounding\' phantom_id '_imageQuality'])

min_specSNR = 5;
min_scatSNR = 1.25;
min_CNR = 0.9;

enough_specSNR = specSNR > min_specSNR;
enough_scatSNR = scatSNR > min_scatSNR;
enough_scatSNR(1,1) = 0;
enough_CNR = CNR > min_CNR;

image_good_enough = enough_specSNR + enough_scatSNR + enough_CNR;

plot_fields = {enough_specSNR, enough_scatSNR, enough_CNR, image_good_enough};
plot_titles = {['specSNR > ' num2str(min_specSNR)], ...
               ['scatSNR > ' num2str(min_scatSNR)], ...
               ['CNR > ' num2str(min_CNR)], ...
               'image good enough' };
img_files = {['specSNR_gt_' num2str(min_specSNR)], ...
             ['scatSNR_gt_' num2str(min_scatSNR)], ...
             ['CNR_gt_' num2str(min_CNR)], ...
             'image_good_enough' };

display_and_save_multiple_fields_vs_bw_cf(plot_fields, plot_titles, img_files)


%% compounding
% uses image_good_enough binary field to choose images for compounding

tic

Nx = 172; Ny = 171; Nt = 961;
compound_image = zeros([Nx,Ny,Nt]);
num_images_compounded = 0;

bw_index = 1;
for bandwidth = bandwidths

    cf_index = 1;
    for centre_freq = centre_freqs

        if image_good_enough(bw_index,cf_index)
            reflection_image = load_image(phantom_id, centre_freq, bandwidth);
            compound_image = compound_image + reflection_image;
            num_images_compounded = num_images_compounded + 1;
        end

        cf_index = cf_index + 1;
    end
    bw_index = bw_index + 1;
end

compound_image = compound_image / num_images_compounded;

disp(['  completed in ' scaleTime(toc)]);

volume_data = reshape(compound_image,Nx,Ny,Nt);
volume_spacing = [kgrid.dx, kgrid.dy, params.dt*c0];

% file_path = ['recon_data\' phantom_id '_compound.mat'];
file_path = ['recon_data\' phantom_id '_compound_bw10.mat'];
save(file_path,'volume_data','volume_spacing','-v7.3')

sliceViewer


%% assess quality of unfiltered, compound and compound with bw10 only

% file_path = ['recon_data\' phantom_id '.mat'];
% file_path = ['recon_data\' phantom_id '_compound.mat'];
file_path = ['recon_data\' phantom_id '_compound_bw10.mat'];
image_data = load(file_path);
compound_image = image_data.volume_data;

% compound_image = compound_image(:,:,80:1040); % trim unfiltered image to same depth

meanIP = get_mean_intensity_projection(compound_image);

signal_tube = get_peak_signal_of_tube(meanIP);
[scatter_water_mean, scatter_water_std] = get_scattering_distr_in_water(meanIP);
[scatter_atmm_mean, scatter_atmm_std] = get_scattering_distr_in_atmm(meanIP);

specSNR = signal_tube       / scatter_water_mean;
scatSNR = scatter_atmm_mean / scatter_water_mean;
CNR = (scatter_atmm_mean - scatter_water_mean) / (scatter_atmm_std + scatter_water_std);


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

function [scatter_water_mean, scatter_water_std] = get_scattering_distr_in_water(meanIP) % hard numbers !
% limits customised to third tube

    ROI = meanIP(61:65,385:425);
    
    scatter_water_mean = mean(ROI(:));
    scatter_water_std  = std(ROI(:));

end

function [scatter_atmm_mean, scatter_atmm_std] = get_scattering_distr_in_atmm(meanIP) % hard numbers !
% limits customised to third tube

    ROI_left  = meanIP(51:55,385:425);
    ROI_right = meanIP(72:76,385:425);
    ROI_front = meanIP(61:65,300:340);
    ROI_back  = meanIP(61:65,465:505);
    
    ROI = [ROI_left ROI_right ROI_front ROI_back];

    scatter_atmm_mean = mean(ROI(:));
    scatter_atmm_std  = std(ROI(:));

end

function display_and_save_multiple_fields_vs_bw_cf(plot_fields, plot_titles, img_files)

    idx = 1;
    for plot_field = plot_fields

        figure
        imagesc(plot_field{1})
            title(plot_titles{idx})
            xlabel('centre frequency [MHz]')
            ylabel('bandwidth [MHz]')
            yticks(1:13)
            yticklabels({'1','2','3','4','6','8','10','15','20','25','30','35','40'})
            run('formatFigures.m')

        file_img = ['..\figures\_Matlab figs\freqCompounding\' img_files{idx} ];
        saveas(gcf,[file_img '.fig'])
        saveas(gcf,[file_img '.jpg'])

        idx = idx + 1;
    end

end


