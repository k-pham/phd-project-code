% USimagingPhantoms.m:      parameters and file locations
% USimagingRecon.m:         run reconstruction, save .mat files, plot
% USimagingAnalysis.m:      image analysis


%% get phantom parameters and SGL file locations from USimagingPhantoms.m

clear all %#ok<CLALL>
run('USimagingPhantoms.m')


%% run multiple reconstructions in loop (start)

% play with c0 & t0 correction:
% for samples_t0_correct = -9
% for c0 = 1477:1:1481

% multiple file names:
scanIDs = 1:length(file_names);
for scanID = scanIDs(1:end)
    
    file_name = file_names{scanID};
    trigger_delay = trigger_delays{scanID};

%% load & view SGL data

display(['Viewing: ' file_name])

[sensor_data, params] = loadSGL([file_dir file_name]);
% sensor_data = sensor_data(:,:,1:1400);
params.trigger_delay        = trigger_delay;
params.Nt_zero_pad_source   = samples_cut_off;
params.Nt_t0_correct        = samples_t0_correct;
params.file_data            = file_name;

% dim = 2;
% sensor_data = squeeze(sensor_data(:,47,:));

figure
set(gcf,'Position',[100,50,600,800])
switch dim
    case 2
        imagesc(sensor_data')
    case 3
        imagesc(squeeze(sensor_data(:,50,:))')
end
    colormap(gray)
    colorbar
    title(strtok(file_name,'@'),'Interpreter','None')
    xlabel('x axis [dx]')
    ylabel('time [dt]')
    % ylim([1,50])

% pause


%% run reconstruction

display(['Reconstructing: ' file_name])

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
                                'SaveImageToFile', true ...
                            );
end


%%

% sliceViewer


%% plot reconstructed image

global kgrid t_array Nt %#ok<TLEV>

figure%(1); clf(1)
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
    cmap = flipud(cmap);    % flip colormap to make black = signal
    colormap(cmap);
    colorbar
    %hold on
    %plot([0,0],[-10,10],'b--')
    %legend('sensor')
    %hold off
%     ylim([40,140])
%     xlim([1250,1450])
%     ylim([0.1,0.8])
%     xlim([0.5,2.5])


%% find image peaks and FWHM for resolution measurements

threshold = 1000;
peaksInfo = imagePeakFinder(reflection_image, c0, threshold);

% concatenate peaksInfo array for all line scans
if(~exist('peaksInfoAll','var'))
    peaksInfoAll = peaksInfo;
else
    peaksInfoAll = cat(2,peaksInfoAll,peaksInfo);
end

%% run multiple reconstructions in loop (end)

% pause

end     % of scanID / file_name loop

% end     % of c0 loop
% end     % of samples_t0_correct loop



%% plot reconstructed image

% figure(9); clf(9)
% Nx = size(reflection_image,1);
% imagesc([180,220]*dt*c0, [180,220]*dx*1e3-Nx/2*dx*1e3, reflection_image(180:220,180:220)) % 
% % 1st index (x) = row index = y axis
%     title('reconstructed image')
%     xlabel('z [mm]')
%     ylabel('x [mm]')
%     %axis images
%     colormap(cmap);
%     colorbar


%% plot MIP (line) of reconstructed image

% reflection_image_MIP = squeeze(max(reflection_image(:,:),[],2));      % p_xz to p.max(z)_x
% 
% figure%(10); clf(10)
% plot( reflection_image_MIP)
%     title('MIP of reconstructed image')
%     xlabel('x [mm]')
%     ylabel('signal amplitude [a.u.]')
% %     xlim([-11.4,11.4])


%% line profile through reconstructed image at given depth for LSF

% figure(11); clf(11)

% slice = 210; % 192, 198, 210
% depth = slice * dt * c0 *1e3 ; % in mm

% plot(kgrid.x_vec*1e3, squeeze(reflection_image(:,slice))) %
%     title(['profile through reconstructed 2D image at depth z = ' num2str(depth) ' mm' ])
%     xlabel('x [mm]')



