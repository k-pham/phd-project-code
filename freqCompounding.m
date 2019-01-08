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


%% generate lots of different reconstructions with varying freq filters
% and save image in .mat and meanIP in .jpg/.fig

bandwidths   = [1:1:3,4:2:10,15:5:40] *1e6;
centre_freqs = (1:1:35) *1e6;

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


%% compounding

centre_freq = [5:5:35]*1e6;
bandwidth = 2e6;

compound_image = zeros([Nx,Ny,Nz]);

for f = centre_freq

    file_path = ['recon_data\' phantom_id ...
         '_f' num2str(f/1e6) ...
         '_bw' num2str(bandwidth/1e6) ...
         '.mat'];
    data = load(file_path);
    
    compound_image = compound_image + data.volume_data;

end

compound_image = compound_image / length(centre_freq);

volume_data = reshape(compound_image,Nx,Ny,Nz);
volume_spacing = [kgrid.dx, kgrid.dy, params.dt*c0];

file_path = ['recon_data\' phantom_id '_compound.mat'];
save(file_path,'volume_data','volume_spacing','-v7.3')

sliceViewer


%% local functions

function display_and_save_meanIP(reflection_image,c0,centre_freq,bandwidth,phantom_id)

    global kgrid t_array

    meanIP = squeeze(mean(reflection_image(:,20:40,:),2));
    
    figure
    imagesc(kgrid.x_vec*1e3,t_array*c0*1e3,meanIP')
        axis image
        colormap(gray)
%         caxis([0 ])
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









