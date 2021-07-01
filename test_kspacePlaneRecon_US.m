
file_dir = '..\data\imagingUS\';

file_data = '180626\polymerLeaf2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_51s28m20h_26-06-18_avg1_2D_raw.SGL';
trigger_delay = 4e-6;

dim = 3;

samples_cut_off = 0;
samples_t0_correct = -6;
c0 = 1484;

disp(['Loading: ' file_data])

[sensor_data, params] = loadSGL([file_dir file_data]);

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_data;

disp(['Reconstructing: ' file_data])

[reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
                            'ZeroPad', 10, ...
                            'Upsample', false, ...
                            'Apodise', false, ...
                            'FreqBandFilter', {}, ...
                            'FreqLowFilter', {}, ...
                            'TimeGainCompensate', {}, ...
                            'EnvelopeDetect', true, ...
                            'LogCompress', 0, ...
                            'SaveImageToFile', false ...
                        );

%% save control image data set
% disp('Saving data to .mat ...'),
% tic
% 
% global Nx Ny kgrid
% Nz = size(reflection_image,3);
% volume_data = reshape(reflection_image,Nx,Ny,Nz);
% volume_spacing = [kgrid.dx, kgrid.dy, params.dt*c0];
% 
% phantom_id = strtok(file_data,'@'); % parse string up to specified delimiter
% phantom_id = phantom_id(8:end);     % remove date folder from string
% 
% file_image = ['recon_data\' phantom_id '_control.mat'];
% 
% save(file_image,'volume_data','volume_spacing','-v7.3')
% 
% disp(['  completed in ' scaleTime(toc)]);
% 
% % sliceViewer

%% load control image data set

disp('Loading control image')

file_control    = 'D:\PROJECT\code\recon_data\polymerLeaf2_BK31[CNT]_control.mat';
control         = load(file_control);
control_image   = control.volume_data;
control_spacing = control.volume_spacing;

%% compare image with control

disp('Comparing to control image')

difference = control_image - reflection_image;
% disp(max(difference,[],'all'))
% disp(min(difference,[],'all'))
if any(difference(:)) == false
    disp('Test passed :)')
else
    disp('Test failed :(')
end

