%% recon params

c = 1484; % [m/s]


%% load data

file_dir  = 'D:\PROJECT\data\angleComp\';

% load data for wire 2D scans on bench-top scanner
% --------------------------------------
file_data_list = {
'210827\wire_1D_angledCNT_8@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_20s20m19h_27-08-21_avg1_1D_raw.SGL'
'210827\wire_1D_angledCNT_6@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_01s19m19h_27-08-21_avg1_1D_raw.SGL'
'210827\wire_1D_angledCNT_4@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_37s17m19h_27-08-21_avg1_1D_raw.SGL'
'210827\wire_1D_angledCNT_2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_55s15m19h_27-08-21_avg1_1D_raw.SGL'
'210827\wire_1D_angledCNT_0@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_01s14m19h_27-08-21_avg1_1D_raw.SGL'
'210827\wire_1D_angledCNT_-2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_59s21m19h_27-08-21_avg1_1D_raw.SGL'
'210827\wire_1D_angledCNT_-4@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_35s23m19h_27-08-21_avg1_1D_raw.SGL'
'210827\wire_1D_angledCNT_-6@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_48s24m19h_27-08-21_avg1_1D_raw.SGL'
'210827\wire_1D_angledCNT_-8@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_17s26m19h_27-08-21_avg1_1D_raw.SGL'
};
file_data_list = reshape(file_data_list,[1 length(file_data_list)]);

for file_data_cell = file_data_list

    file_data = file_data_cell{1};
    id = strtok(file_data,'@');
	id = str2double(id(26:end));
    
    disp(['angle id: ' num2str(id)])
% --------------------------------------


%%
[sensor_data, sensor_params] = loadSGL([file_dir file_data]);

y_line = round(sensor_params.Ny/2);

fig_data1 = figure(1);
imagesc(sensor_data')
    colormap(gray)
    colorbar
    title(strtok(file_data,'@'),'Interpreter','None')
    xlabel('y axis [dx]')
    ylabel('time [dt]')
    drawnow

fig_data2 = figure(2);
plot(squeeze(sensor_data(y_line,:)))
    xlabel('time [dt]')
    ylabel('y axis [dx]')


%% extract TOA from timeseries data

[~,TOA_y_idx] = max(sensor_data,[],2);     % [time index]
TOA_y = TOA_y_idx * sensor_params.dt;      % [s]

kgrid = kWaveGrid(sensor_params.Nx,sensor_params.dx,sensor_params.Ny,sensor_params.dy);
y_vec = kgrid.y_vec;


%% initial fit of plane wave to source

% regression in 1d TOA(x) = p1 + p2*x
Model = [ones(length(y_vec),1) y_vec];
params = Model \ TOA_y;

% evaluate fit over whole scanning area
TOA_fit = params(1) + params(2) * y_vec;

% get deviation in time index of arrival from fit within mask
deviation = TOA_y - TOA_fit;


%% exclude points with TOA too far away from fit

num_dt_threshold = 5;
pointsToKeep = deviation < num_dt_threshold*sensor_params.dt;
TOA_y        = TOA_y(pointsToKeep);
y_vec        = y_vec(pointsToKeep);


%% redo fit of plane wave to source without excluded points

% regression in 1d TOA(x) = p1 + p2*x
Model = [ones(length(y_vec),1) y_vec];
params = Model \ TOA_y;

% evaluate fit over whole scanning area
TOA_fit = params(1) + params(2) * y_vec;

% get deviation in time index of arrival from fit within mask
deviation = TOA_y - TOA_fit;


%% plot fits and deviation

% plot deviation of TOA from TOA_fit in 2d colour scatter
figure(3)
plot(y_vec*1e3,deviation*1e9)
% hold on
    title('time of arrival: deviation from fit within ROI [s]')
    xlabel('y axis / mm')
    ylabel('deviation from fit [ns]')
    
% plot 3d scatter of TOA for whole scanning area with fit from masked data
figure(4)
plot(y_vec*1e3,TOA_y*1e6,'.')
hold on
plot(y_vec*1e3,TOA_fit*1e6)
    title('Time of arrival for each interrogation point with fit')
    xlabel('y axis / mm')
    ylabel('time of arrival / \mus')
    drawnow
% hold off


%% calculate steering angle components ax and ay
% TOA(x) = p1 + p2*x
% c*p2 = sin(a)

a = asin(c*params(2));

disp(['steering angle of plane wave is: ' num2str(rad2deg(a)) ' deg'])
disp('=====================================')


%%
% pause
% end % angles
% pause

%% t0 correction

% indeces for central timeseries
nyc = round(sensor_params.Ny/2);

% get TOA(source) in central timeseries
TOA_source = TOA_y_idx(nyc);

%% set t0(excitation) in dt
for t0_excitation = 8 % 5:1:11 % -15:5:20

disp(['t0(excitation) = ' num2str(t0_excitation) '*dt'])

% reload sensor_data and _params to avoid reapplying t0 correction
[sensor_data, sensor_params] = loadSGL([file_dir file_data]);


%% remove acoustic source from signal

source_pads = 2000;
sensor_data = cat(2, zeros(kgrid.Ny,source_pads), sensor_data(:,source_pads+1:end));

% check that padding did not change length of sensor_data
assert(size(sensor_data,2) == sensor_params.Nt)


%% zero padding source offset to sensor plane
% assumes that correction involves zero padding rather than trimming

% set t0(source) in dt
t0_source = t0_excitation - (TOA_source - t0_excitation);

% number of pads
pads = -t0_source;

% check that correction requires zero-padding and not trimming
assert(pads>0)

% apply t0 correction
sensor_data_padded = cat(2, zeros(sensor_params.Ny,pads), sensor_data);
sensor_params.Nt   = sensor_params.Nt + pads;


%% reconstruction

reflection_image = kspaceLineRecon_US_steered(sensor_data_padded',sensor_params.dy,sensor_params.dt,c,a);
% reflection_image = kspaceLineRecon_US(sensor_data_padded',sensor_params.dy,sensor_params.dt,c);
                                                        % input p_tx, output p_zx
reflection_image = permute(reflection_image,[2 1]);     % reorder p_zx to p_xz
reflection_image = reflection_image(:,1:round(size(reflection_image,2)/2));  % trim half of z if using doubled depth kgrid


%% envelope detection

disp('Envelope detecting ...')
tic
reflection_image_env = envelopeDetection(reflection_image);
% reflection_image_env = transpose(envelopeDetection(reflection_image'));
assert(isequal( size(reflection_image_env), size(reflection_image) ))
disp(['  completed in ' scaleTime(toc)]);


%% plot image

t_array = sensor_params.dt * (1:1:sensor_params.Nt);
z_vec   = t_array * c / 2;

figure
imagesc(kgrid.y_vec*1e3,z_vec*1e3,reflection_image_env')
title(['angle ' num2str(id) ' ' num2str(rad2deg(a)) ', t0 = ' num2str(t0_excitation) '*dt'])
xlabel('y axis / mm')
ylabel('depth / mm')
xlim([-1.1,-0.9])
ylim([10.1,10.5])
drawnow


%% peak analysis of wire

wire_peak = max(reflection_image_env,[],'all');
disp(['MAX OF WIRE = ' num2str(wire_peak)])


%%
disp('=====================================')
% pause
end % t0
end % angles



