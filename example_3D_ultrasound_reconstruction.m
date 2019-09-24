%%

file_dir = '..\data\imagingUS\';
file_name = '180626\optifibreKnot_angled4_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_10s39m17h_26-06-18_avg1_2D_raw.SGL';

dim = 3;

c0 = 1484;

trigger_delay = 5e-6;
samples_cut_off = 0;
samples_t0_correct = -6;


%% load SGL data & add parameters to params

[sensor_data, params] = loadSGL([file_dir file_name]);

params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;


%% view SGL data

fig_data = figure;
set(gcf,'Position',[100,50,600,800])
switch dim
    case 2
        imagesc(sensor_data')
    case 3
        imagesc(squeeze(sensor_data(40,:,:))')
end
    colormap(gray)
    colorbar
    title(strtok(file_name,'@'),'Interpreter','None')
    xlabel('x axis [dx]')
    ylabel('time [dt]')
    drawnow


%% run reconstruction

if dim == 2
    [reflection_image] = reconstruct2dUSimage(sensor_data, params, c0);
elseif dim == 3
    [reflection_image] = reconstruct3dUSimage(sensor_data, params, c0, ...
                                'ZeroPad', 10, ...
                                'Upsample', false, ...
                                'Apodise', false, ...
                                'FreqBandFilter', {10e6,100e6}, ...
                                'TimeGainCompensate', {}, ...
                                'EnvelopeDetect', true, ...
                                'LogCompress', 0, ...
                                'SaveImageToFile', false ...
                            );
end


%%

% sliceViewer


%% plot reconstructed image

global kgrid t_array Nt

fig_image = figure;
set(gcf,'Position',[100,100,800,450])
switch dim
    case 2
        % imagesc(reflection_image(900:1020,160:280)')
%         imagesc(reflection_image(:,1:samples_total/2)')
        imagesc(kgrid.x_vec*1e3, t_array*c0*1e3, reflection_image(:,1:Nt)') % omit factor 1/2 in dz because of doubled depth bug
            % 1st index (x) = row index = y axis -- transposed -> x axis
            xlabel('x [mm]')
            ylabel('z [mm]')
    case 3
        imagesc(kgrid.x_vec*1e3, t_array*c0*1e3, squeeze(reflection_image(:,75,1:Nt))') % omit factor 1/2 in dz because of doubled depth bug
            xlabel('x [mm]')
            ylabel('z [mm]')
end
    title(['reconstructed image with c0 = ' num2str(c0) ', t0 correction = ' num2str(params.Nt_t0_correct)])
    %axis image
    %caxis([0,200])
    cmap = colormap(gray);
%     cmap = flipud(cmap);    % flip colormap to make black = signal
    colormap(cmap);
    colorbar
    drawnow
    %hold on
    %plot([0,0],[-10,10],'b--')
    %legend('sensor')
    %hold off
%     ylim([40,140])
%     xlim([1250,1450])
%     ylim([0.1,0.8])
%     xlim([0.5,2.5])




