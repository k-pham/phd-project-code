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
% phantom_id = 'porkBelly3_BK31[CNT]';
% file_name = '180828\porkBelly3_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_26s20m19h_28-08-18_avg1_2D_raw.SGL';
% trigger_delay = 0;
% samples_cut_off = 10;
% samples_t0_correct = -4;
% c0 = 1460;

% 190114 lymph node (L3)
% phantom_id = 'lymphNode2_BK31[CNT]';
% file_name = '190114/lymphNode2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_13s51m16h_14-01-19_avg1_2D_raw.SGL';
% trigger_delay = 0;
% samples_cut_off = 50;
% samples_t0_correct = -6;
% c0 = 1520;


%% frequency bands to compound

bandwidths   = 10e6;
centre_freqs = (2:1:15)*1e6;


%% load and prepare data

[sensor_data, params] = loadSGL([file_dir file_name]);

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;

global kgrid t_array %#ok<NUSED>


%% set up compounding

Nx = 344; Ny = 342; Nt = 1497;
compound_image = zeros([Nx,Ny,Nt]);
num_images_compounded = 0;


%% generate lots of different reconstructions with varying freq filters
% and save image in .mat and meanIP in .jpg/.fig

for bandwidth = bandwidths
    for centre_freq = centre_freqs

        disp(['CENTRE FREQ ' num2str(centre_freq/1e6) ', BANDWIDTH ' num2str(bandwidth/1e6)])
        
        [reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
                                    'ZeroPad', 10, ...
                                    'Upsample', true, ...
	                                'Apodise', false, ...
                                    'FreqBandFilter', {centre_freq, bandwidth}, ...
                                    'FreqLowFilter', {}, ...
                                    'TimeGainCompensate',{'Exponential',50}, ...
                                    'EnvelopeDetect', true, ...
                                    'LogCompress', 0, ...
                                    'SaveImageToFile', true ...
                                );
                            
        display_and_save_meanIP(reflection_image, c0, centre_freq, bandwidth, phantom_id)
                
        compound_image = compound_image + reflection_image;
        num_images_compounded = num_images_compounded + 1;
        
    end
end

compound_image = compound_image / num_images_compounded;

display_and_save_meanIP(compound_image,c0,2.15e6,bandwidth,phantom_id)


%% local functions

function meanIP = get_mean_intensity_projection(reflection_image) % hard numbers !
% limits customised to atmm_orgasol1_BK31[CNT]

    global t_array
    
    meanIP = squeeze(mean(reflection_image(:,40:80,80:1040),2));
    t_array = t_array(80:1040);

end

function display_and_save_meanIP(reflection_image,c0,centre_freq,bandwidth,phantom_id,varargin)
% varargin to allow for optional specification of colour axis

    global kgrid t_array

    meanIP = get_mean_intensity_projection(reflection_image);
    
    figure
    imagesc(kgrid.x_vec*1e3,t_array*c0*1e3,meanIP')
        axis image
        colormap(gray)
        if ~isempty(varargin)
            caxis(varargin{1})
        end
        brighten(0.5)
        title(['z-x MeanIP 2 mm slice: f = ' num2str(centre_freq/1e6) ', bw = ' num2str(bandwidth/1e6)])
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

