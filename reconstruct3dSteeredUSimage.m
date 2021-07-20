%% load data

file_dir  = 'D:\PROJECT\data\angleComp\210630\';
file_data = 'polymerLeaf_angledCNT_0@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_26s35m20h_30-06-21_avg1_2D_raw.SGL';

[sensor_data, params] = loadSGL([file_dir file_data]);

fig_data = figure;
set(gcf,'Position',[100,50,600,800])
imagesc(squeeze(sensor_data(70,:,:))')
    colormap(gray)
    colorbar
    title(strtok(file_data,'@'),'Interpreter','None')
    xlabel('x axis [dx]')
    ylabel('time [dt]')
    drawnow
fig_data2 = figure;
set(gcf,'Position',[500,50,800,600])
plot(squeeze(sensor_data(70,70,:)))
    xlabel('time [dt]')
    ylabel('x axis [dx]')


%% fit plane wave to source



%% calculate steering angle (a), rotation angle (b) & time offset (t0)



%% t0 correction



%% zero padding delay


