%% recon params

c = 1484; % [m/s]


%% load data

file_dir  = 'D:\PROJECT\data\angleComp\';
file_data = '210630\polymerLeaf_angledCNT_0@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_26s35m20h_30-06-21_avg1_2D_raw.SGL';
% file_data = '210827\wire_angledCNT_0@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_26s04m19h_27-08-21_avg1_2D_raw.SGL';

% load data for wire 3D scans on 16-beam
% --------------------------------------
% file_data_list = {
% '210818\Full_scan1_wire27_CNT_angle0@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_41s11m16h_18-08-21_avg1_2D_raw.SGL'
% '210818\Full_scan1_wire27_CNT_angle2@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_28s14m16h_18-08-21_avg1_2D_raw.SGL'
% '210818\Full_scan1_wire27_CNT_angle4@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_53s17m16h_18-08-21_avg1_2D_raw.SGL'
% '210818\Full_scan1_wire27_CNT_angle6@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_51s20m16h_18-08-21_avg1_2D_raw.SGL'
% '210818\Full_scan1_wire27_CNT_angle8@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_37s24m16h_18-08-21_avg1_2D_raw.SGL'
% '210818\Full_scan1_wire27_CNT_angle-2@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_58s26m16h_18-08-21_avg1_2D_raw.SGL'
% '210818\Full_scan1_wire27_CNT_angle-4@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_08s29m16h_18-08-21_avg1_2D_raw.SGL'
% '210818\Full_scan1_wire27_CNT_angle-6@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_40s32m16h_18-08-21_avg1_2D_raw.SGL'
% '210818\Full_scan1_wire27_CNT_angle-8@850nm_t0[0]_dx[109µm]_dy[108µm]_dt[17ns]_37s35m16h_18-08-21_avg1_2D_raw.SGL'
% };
% --------------------------------------

% load data for wire 3D scans on bench-top scanner
% --------------------------------------
file_data_list = {
'210906\wire_angledCNT_8@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_05s39m17h_06-09-21_avg1_2D_raw.SGL'
'210906\wire_angledCNT_4@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_35s18m17h_06-09-21_avg1_2D_raw.SGL'
'210906\wire_angledCNT_0@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_04s47m16h_06-09-21_avg1_2D_raw.SGL'
'210906\wire_angledCNT_-4@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_54s59m17h_06-09-21_avg1_2D_raw.SGL'
'210906\wire_angledCNT_-8@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_28s20m18h_06-09-21_avg1_2D_raw.SGL'
};
% --------------------------------------

file_data_list = reshape(file_data_list,[1 length(file_data_list)]);

for file_data_cell = file_data_list(1)

    file_data = file_data_cell{1};
    id = strtok(file_data,'@');
	id = str2double(id(23:end)); % 35
    
    disp(['angle id: ' num2str(id)])


%%
[sensor_data, sensor_params] = loadSGL([file_dir file_data]);

x_slice = round(sensor_params.Nx/2);
y_slice = round(sensor_params.Ny/2);

fig_data1 = figure(1);
% set(gcf,'Position',[100,50,600,800])
imagesc(squeeze(sensor_data(:,y_slice,:))')
    colormap(gray)
    colorbar
    title(strtok(file_data,'@'),'Interpreter','None')
    xlabel('x axis [dx]')
    ylabel('time [dt]')
    drawnow
fig_data2 = figure(2);
% set(gcf,'Position',[100,50,600,800])
imagesc(squeeze(sensor_data(x_slice,:,:))')
    colormap(gray)
    colorbar
    title(strtok(file_data,'@'),'Interpreter','None')
    xlabel('y axis [dx]')
    ylabel('time [dt]')
    drawnow
fig_data3 = figure(3);
% set(gcf,'Position',[500,50,800,600])
plot(squeeze(sensor_data(x_slice,y_slice,:)))
    xlabel('time [dt]')
    ylabel('x axis [dx]')


%% extract TOA from timeseries data

% note peaks will be saturated
% SOLUTION: pick first index at max val as TOA
% ALTERNATIVE: pick first zero crossing (but need to subtract DC offset)
% IN PRACTICS: may need to keep this a variable relative to first max

% get PMAX and TOA in 2d grid and vector (if max element occurs more than 
% once, output = index to 1st occurrence)
[~,TOA_xy_idx] = max(sensor_data,[],3);     % [time index]

% convert TOA from index to realtime
TOA_xy = TOA_xy_idx * sensor_params.dt;     % [s]


%% reshape to vectors

kgrid = kWaveGrid(sensor_params.Nx,sensor_params.dx,sensor_params.Ny,sensor_params.dy);

x = kgrid.x_vec;
y = kgrid.y_vec;
[X,Y] = meshgrid(x,y); X=X'; Y=Y';

Xv = reshape(X,[kgrid.Nx*kgrid.Ny 1]);
Yv = reshape(Y,[kgrid.Nx*kgrid.Ny 1]);

TOA_vec  = reshape(TOA_xy, [kgrid.Nx*kgrid.Ny 1]);


%% maybe apply mask for central area to avoid outliers at edge ruining fit

apply_mask = false;

switch apply_mask
    case true
        
        % check if there is a variable mask, and if not
        if ~exist('mask','var')
            
            % check if there is a saved mask
            file_mask = [file_dir file_data(1:7) 'mask.mat'];
            
            % if yes, load mask
            if exist(file_mask,'file')
                load(file_mask,'mask')
                
            % if not, make mask & save
            else
                canvas = ones(kgrid.Nx,kgrid.Ny);
                figure
                imagesc(x,y,canvas)
                mask = roipoly;
                save(file_mask,'mask')
            end
            
        % if yes, use existing mask
        else
            disp('Using existing mask')
        end
        
        % apply mask
        TOA_vec_mask = TOA_xy(mask);
        Xv_mask = Xv(mask);
        Yv_mask = Yv(mask);
        
    % if no mask applied do this to make into vectors
    case false
        TOA_vec_mask = TOA_vec;
        Xv_mask = Xv;
        Yv_mask = Yv;
end


%% initial fit of plane wave to source

% regression in 2d TOA(x,y) = p1 + p2*x + p3*y using masked data points
Model = [ones(length(Xv_mask),1) Xv_mask Yv_mask];
params = Model \ TOA_vec_mask;

% evaluate fit over whole scanning area
TOA_fit_vec = params(1) + params(2) * Xv + params(3) * Yv;
%TOA_fit_xy = reshape(TOA_fit_vec,[kgrid.Nx kgrid.Ny]);

% get deviation in time index of arrival from fit within mask
TOA_fit_vec_mask = TOA_fit_vec;%(mask);
deviation_mask = TOA_vec_mask - TOA_fit_vec_mask;


%% exclude points with TOA too far away from fit

num_dt_threshold = 5;
pointsToKeep = deviation_mask < num_dt_threshold*sensor_params.dt;
TOA_vec_mask = TOA_vec_mask(pointsToKeep);
Xv_mask      = Xv_mask(pointsToKeep);
Yv_mask      = Yv_mask(pointsToKeep);


%% redo fit of plane wave to source without excluded points

% regression in 2d TOA(x,y) = p1 + p2*x + p3*y using masked data points
Model = [ones(length(Xv_mask),1) Xv_mask Yv_mask];
params = Model \ TOA_vec_mask;

% evaluate fit over whole scanning area
TOA_fit_vec = params(1) + params(2) * Xv + params(3) * Yv;
TOA_fit_xy = reshape(TOA_fit_vec,[kgrid.Nx kgrid.Ny]);

% get deviation in time index of arrival from fit within mask
TOA_fit_vec_mask = TOA_fit_vec;%(mask);
TOA_fit_vec_mask = TOA_fit_vec_mask(pointsToKeep);
deviation_mask = TOA_vec_mask - TOA_fit_vec_mask;

% plot deviation of TOA from TOA_fit in 2d colour scatter
figure(4)
scatter(Xv_mask,Yv_mask,[],deviation_mask,'o')
% hold on
    colorbar
    %caxis([-3,3])
    title('time of arrival: deviation from fit within ROI [s]')
    xlabel('x')
    ylabel('y')
    %xlim([30,120])
    %ylim([20,120])
    axis equal
    
% plot 3d scatter of TOA for whole scanning area with fit from masked data
figure(5)
scatter3(Xv_mask,Yv_mask,TOA_vec_mask,'.')
hold on
mesh(X,Y,TOA_fit_xy)
    title('Time of arrival for each interrogation point with fit')
    xlabel('x / m')
    ylabel('y / m')
    zlabel('time of arrival / s')
    drawnow
%     zlabel('time index of arrival')
%     xlim([0,kgrid.Nx])
%     ylim([0,kgrid.Ny])
%     zlim([160,240])
% zlim([0,0.002])


%% calculate steering angle components ax and ay
% TOA(x,y) = p1 + p2*x + p3*y
% c*p2 = sin(ax)
% c*p3 = sin(ay)

ax = asin(c*params(2));
ay = asin(c*params(3));

disp(['angle of plane wave in x direction is: ' num2str(ax/pi*180) ' deg'])
disp(['angle of plane wave in y direction is: ' num2str(ay/pi*180) ' deg'])


%% calculate steering angle (a), rotation angle (b)
% tan(a) = sqrt( tan^2(ax) + tan^2(ay) )
% tan(b) = tan(ax) / tan(ay)

a = atan( sqrt( tan(ax).^2 + tan(ay).^2 ) );
b = atan( tan(ay) ./ tan(ax) );

disp(['steering angle of plane wave is: ' num2str(a/pi*180) ' deg'])
disp(['rotation angle of plane wave is: ' num2str(b/pi*180) ' deg'])


%%
% pause
% end


%% t0 correction

% indeces for central timeseries
nxc = round(sensor_params.Nx/2);
nyc = round(sensor_params.Ny/2);

% get TOA(source) in central timeseries
TOA_source = TOA_xy_idx(nxc,nyc);

%% set t0(excitation) in dt
for t0_excitation = 0 %-15:5:20

disp(['t0(excitation) = ' num2str(t0_excitation) '*dt'])

% reload sensor_data and _params to avoid reapplying t0 correction
[sensor_data, sensor_params] = loadSGL([file_dir file_data]);


%% remove acoustic source from signal

source_pads = 2000;
sensor_data = cat(3, zeros(kgrid.Nx,kgrid.Ny,source_pads), sensor_data(:,:,source_pads+1:end));

% check that padding did not change length of sensor_data
assert(size(sensor_data,3) == sensor_params.Nt)


%% zero padding source offset to sensor plane
% assumes that correction involves zero padding rather than trimming

% set t0(source) in dt
t0_source = t0_excitation - (TOA_source - t0_excitation);

% number of pads
pads = -t0_source;

% check that correction requires zero-padding and not trimming
assert(pads>0)

% apply t0 correction
sensor_data_padded = cat(3, zeros(sensor_params.Nx,sensor_params.Ny,pads), sensor_data);
sensor_params.Nt   = sensor_params.Nt + pads;


%% reconstruction

sensor_data_padded = permute(sensor_data_padded,[3 1 2]);   % reorder p_xyt to p_txy

reflection_image = kspacePlaneRecon_US_steered(sensor_data_padded,sensor_params.dx,sensor_params.dy,sensor_params.dt,c,a,b);
% reflection_image = kspacePlaneRecon_US(sensor_data_padded,sensor_params.dx,sensor_params.dy,sensor_params.dt,c);
                                                            % output as p_zxy
reflection_image = permute(reflection_image,[2 3 1]);       % reorder p_zxy to p_xyz


%% envelope detection

disp('Envelope detecting ...')
tic

reflection_image_env = zeros(size(reflection_image));
for i = 1:sensor_params.Ny
    reflection_image_env(:,i,:) = envelopeDetection(squeeze(reflection_image(:,i,:)));
end
assert(isequal( size(reflection_image_env), size(reflection_image) ))

disp(['  completed in ' scaleTime(toc)]);


%% save for sliceViewer

volume_data = reshape(reflection_image_env,sensor_params.Nx,sensor_params.Ny,sensor_params.Nt);
volume_spacing = [sensor_params.dx,sensor_params.dy,sensor_params.dt*c];

phantom_id = strtok(file_data,'@'); % parse string up to specified delimiter
phantom_id = phantom_id(8:end);     % remove date folder from string
file_image = ['recon_data\' phantom_id '_t0_' num2str(t0_excitation) 'dt.mat'];

save(file_image,'volume_data','volume_spacing','-v7.3')

% sliceViewer


%% generate MIP yz and xz

z_range = 1:size(reflection_image_env,3);% 1050;

MIP_yz = squeeze(max(reflection_image_env(:,:,z_range),[],1));
MIP_xz = squeeze(max(reflection_image_env(:,:,z_range),[],2));
z_axis = z_range * sensor_params.dt * c;

figure
sgtitle(['t0(excitation) = ' num2str(t0_excitation) '*dt'])

subplot(1,2,1)
imagesc(kgrid.y_vec'*1e3, z_axis*1e3, MIP_yz')
title('MIP yz')
xlabel('y axis / mm')
ylabel('depth / mm')

subplot(1,2,2)
imagesc(kgrid.x_vec'*1e3, z_axis*1e3, MIP_xz')
title('MIP xz')
xlabel('x axis / mm')
ylabel('depth / mm')

drawnow

%%
% pause
end % t0
end % angles


%% save MIPS

% fignums = 8:15;
% t0_excitations = -15:5:20;
% for i=1:8
% fignum = fignums(i);
% t0_excitation = t0_excitations(i);
% figure(fignum)
% file_fig = [file_dir phantom_id '_t0_' num2str(t0_excitation) 'dt'];
% savefig(fignum,[file_fig '.fig'])
% saveas(fignum,[file_fig '.jpg'])
% end

%% maximum intensity projection in xy

% MIP_xy = max(reflection_image_env(:,:,2000:end),[],3);
% 
% figure
% imagesc(MIP_xy)


%% fly-through in z direction

% figure
% for frame = 1:size(reflection_image_env,3)
%     slice = squeeze(reflection_image_env(:,:,frame));
%     imagesc(slice)
%     title(num2str(frame))
%     drawnow
% end


%% fly-through data in time

% figure
% for frame = 1:size(sensor_data_padded,3)
%     slice = squeeze(sensor_data_padded(:,:,frame));
%     imagesc(slice)
%     title(num2str(frame))
%     drawnow
% end

