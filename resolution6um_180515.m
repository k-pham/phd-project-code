%% 180515 resolution6um (carbon fibres in water) move around in space

file_dir = '..\data\USimaging';

dim = 2;
dx = 0.02e-3;
trigger_delay = 0;

samples_cut_off = 100;
samples_t0_correct = 2;
c0 = 1484;

file_names = {
    '180515\resolution6um_BA21[CNT]_z0.5_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_13s15m21h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z0.5_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_04s24m21h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z1_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_20s32m21h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z1_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_35s40m21h_15-05-18_avg10_1D_raw.SGL';
    '180515\resolution6um_BA21[CNT]_z1.5_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_44s52m21h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z1.5_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_39s00m22h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z2_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_33s08m22h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z2_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_04s16m22h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z2.5_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_15s26m22h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z2.5_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_48s33m22h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z3_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_21s41m22h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z3_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_49s48m22h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z4_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_25s58m22h_15-05-18_avg10_1D_raw.SGL'
    '180515\resolution6um_BA21[CNT]_z5_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_27s08m23h_15-05-18_avg10_1D_raw.SGL'
    };

% need to reshape to be able to loop through
file_names = reshape(file_names,[1 length(file_names)]);

for file_name = file_names

    [reflection_image, samples_total, dt, t_array, kgrid] = ...
                reconstruct2dUSimage(file_dir, file_name{1}, dx, trigger_delay, ...
                                    samples_cut_off, samples_t0_correct, c0);

    figure
    set(gcf,'Position',[100,500,800,450])
    % imagesc(reflection_image(:,1:samples_total/2)')
    imagesc(kgrid.x_vec*1e3, t_array*c0/2*1e3, reflection_image(:,1:samples_total/2)')
    % 1st index (x) = row index = y axis
    title(['reconstructed image with c0 = ' num2str(c0)])
    xlabel('x [mm]')
    ylabel('z [mm]')
    %axis image
    %caxis([0,200])
    cmap = colormap(gray);
    cmap = flipud(cmap);    % flip colormap to make black = signal
    colormap(cmap);
    colorbar
    
    drawnow
end






