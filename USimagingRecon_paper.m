%% from USimagingPhantoms.m

file_dir = '..\data\imagingUS\';


% 180626 optiFibreKnot angled 4
%     file_name = '180626\optifibreKnot_angled4_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_10s39m17h_26-06-18_avg1_2D_raw.SGL';
%     trigger_delay = 5e-6;
%     samples_cut_off = 0;
%     samples_t0_correct = -6;
%     c0 = 1484;

% 180626 polymerLeaf curved
%     file_name = '180626\polymerLeaf_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_39s38m19h_26-06-18_avg1_2D_raw.SGL';
%     trigger_delay = 4.5e-6;
%     samples_cut_off = 0;
%     samples_t0_correct = -6;
%     c0 = 1484;

% 180629 gel wax tmm
%     file_name = '180629\gelwaxLayers_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_39s51m16h_29-06-18_avg1_2D_raw.SGL';
%     trigger_delay = 0;
%     samples_cut_off = 10;
%     samples_t0_correct = -6;
%     c0 = 1470;
%         % clinical scanner:
%         file_name_cli = '180629\clinical scanner\29-06-2018_17-41-25 DICOM +2 higher gain\17-38-33.dcm';

% 181204 atmm orgasol 1
%     file_name = '181204/atmm_orgasol1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
%     trigger_delay = 0;
%     samples_cut_off = 50;
%     samples_t0_correct = -6;
%     c0 = 1544;

% 180828 pork belly 3
    file_name = '180828\porkBelly3_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_26s20m19h_28-08-18_avg1_2D_raw.SGL';
    trigger_delay = 0;
    samples_cut_off = 10;
    samples_t0_correct = -4;
    c0 = 1460;

% 190114 lymph node (L3)
%     file_name = '190114/lymphNode2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_13s51m16h_14-01-19_avg1_2D_raw.SGL';
%     trigger_delay = 0;
%     samples_cut_off = 50;
%     samples_t0_correct = -6;
%     c0 = 1520;


%% load sensor data

[sensor_data, params] = loadSGL([file_dir file_name]);
% sensor_data = - sensor_data;            % flip for trolley scanner

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;


%% view sensor data

disp(['Viewing: ' file_name])

fig_data = figure;
set(gcf,'Position',[100,50,600,800])
imagesc(squeeze(sensor_data(40,:,:))')
    colormap(gray)
    colorbar
    title(strtok(file_name,'@'),'Interpreter','None')
    xlabel('x axis [dx]')
    ylabel('time [dt]')
    drawnow

% pause


%% run reconstruction

disp(['Reconstructing: ' file_name])

[reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
                            'ZeroPad', 10, ...
                            'Upsample', true, ...
                            'Apodise', false, ...
                            'FreqBandFilter', {10e6,15e6}, ...
                            'TimeGainCompensate', {}, ...
                            'EnvelopeDetect', true, ...
                            'LogCompress', 0, ...
                            'SaveImageToFile', true ...
                        );

%% post processing

sliceViewer



