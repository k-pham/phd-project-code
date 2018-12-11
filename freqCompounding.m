clear all %#ok<CLALL>
run('USimagingPhantoms.m')

[sensor_data, params] = loadSGL([file_dir file_name]);

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;


centre_freq = 25e6;
bandwidth = 20e6;

[reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
                            'ZeroPad', 10, ...
                            'FreqBandFilter', {centre_freq, bandwidth}, ...
                            'EnvelopeDetect', true ...
                        );


global Nx Ny kgrid

[Nz] = size(reflection_image,3);
volume_data = reshape(reflection_image,Nx,Ny,Nz);
volume_spacing = [kgrid.dx, kgrid.dy, params.dt*c0];

phantom_id = strtok(file_name,'@'); % parse string up to specified delimiter
phantom_id = phantom_id(8:end);     % remove date folder from string
file_path = ['recon_data\' phantom_id ...
             '_f' num2str(centre_freq/1e6) ...
             '_bw' num2str(bandwidth/1e6) ...
             '.mat'];
save(file_path,'volume_data','volume_spacing','-v7.3')

sliceViewer

% title(['z-x MeanIP 1.5 mm slice: f = ' num2str(centre_freq/1e6) ', bw = ' num2str(bandwidth/1e6)])

%% compounding

phantom_id = 'atmm_orgasol1_BK31[CNT]';

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




