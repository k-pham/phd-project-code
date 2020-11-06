%% compound image

compound_image = zeros(2204,3244);


%% USimagingPhantoms

dim = 2;

file_dir = '..\data\imagingUS\191126\';
files = dir(file_dir);
num_files = length(files);

file_name_stem = 'resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D';
file_name_xendings = {'-1','-05','0','05','1'};

num_lines = length(file_name_xendings);
num_files_inline = 13;

% loop through 5 lines in x
for idx_x = 2 %1:num_lines
    
    file_name_filter_x = [file_name_stem, '_x', file_name_xendings{idx_x}];

    file_names_inline = cell(1,num_files_inline);

    % loop through 13 files along z (per line)
    for idx_z = 1:num_files_inline
        file_name_filter_z = [file_name_filter_x, '_z', num2str(idx_z-1) '@'];
        %look through all files
        for idx = 1:num_files
            file = files(idx);
            %with long enough file names
            if(length(file.name) > 7)
                %ending in raw.SGL
                if strcmp(file.name(end-6:end),'raw.SGL')
                    %add to inline if match filter
                    if strcmp(file.name(1:length(file_name_filter_z)),file_name_filter_z)
                        file_names_inline{idx_z} = file.name;
                    end
                end
            end
        end
    end % of loop through 13 files along z (per line)

    trigger_delays = {0, 0, 1e-6, 3e-6, 5e-6, 6e-6, 8e-6, 9e-6, 10e-6, 12e-6, 13e-6, 14e-6, 16e-6};

    samples_cut_off = 50;
%     samples_t0_correct = -14;
%     samples_t0_correct = 2;
    

    %% USimaging Phantom -> Recon

    file_names = file_names_inline;


    %% USimagingRecon

    % save figures to directory path
    dir_figures = 'D:\PROJECT\figures\_Matlab figs\USimaging\191126 resolution27umPlanar BK31[CNT] trolley scrambled fibre centralised parallel phantom\';
    dir_figures = [dir_figures 'x' file_name_xendings{idx_x} ' c0_1489 t0_-14\'];

    % loop through range of c0
    for samples_t0_correct = -14 %-18:-12
    for c0 = 1489 %1484:1494

        % multiple file names:
        scanIDs = 1:length(file_names);
        for scanID = scanIDs(1:end)

            file_name = file_names{scanID};
            trigger_delay = trigger_delays{scanID};

            %% load SGL data

            display(['Loading: ' file_name])

            [sensor_data, params] = loadSGL([file_dir file_name]);
            sensor_data = - sensor_data;            % flip for trolley scanner
            sensor_data = sensor_data(:,1:2500);  % cut-off zero-padding from trolley scanner

            params.trigger_delay        = trigger_delay;
            params.Nt_zero_pad_source   = samples_cut_off;
            params.Nt_t0_correct        = samples_t0_correct;
            params.file_data            = file_name;


            %% view SGL data

            % fig_data = figure;
            % set(gcf,'Position',[100,50,600,800])
            % switch dim
            %     case 2
            %         imagesc(sensor_data')
            %     case 3
            %         imagesc(squeeze(sensor_data(40,:,:))')
            % end
            %     colormap(gray)
            %     colorbar
            %     title(strtok(file_name,'@'),'Interpreter','None')
            %     xlabel('x axis [dx]')
            %     ylabel('time [dt]')
            %     % ylim([1,50])
            %     drawnow
            % 
            % % pause


            %% run reconstruction

            display(['Reconstructing: ' file_name])

            [reflection_image] = reconstruct2dUSimage(sensor_data, params, c0, ...
                                        'Upsample', true, ...
                                        'FreqBandFilter', {}, ...
                                        'FreqCustomFilter', {}, ...
                                        'TimeGainCompensate', {}, ...
                                        'EnvelopeDetect', true, ...
                                        'LogCompress', 0, ...
                                        'SaveImageToFile', false ...
                                    );


            %%

            % sliceViewer


            %% plot line horizontal MIP profile

            % reflection_image_MIP = squeeze(max(reflection_image(:,:),[],2));      % p_xz to p.max(z)_x
            % 
            % fig_profile = figure;
            % plot(reflection_image_MIP)
            %     title('MIP of reconstructed image')
            %     drawnow


            %% plot reconstructed image

            global kgrid %t_array Nt

%             fig_image = figure;
%             set(gcf,'Position',[100,100,800,450])
%             switch dim
%                 case 2
%                     imagesc(kgrid.x_vec, t_array*c0, reflection_image(:,1:Nt)') % omit factor 1/2 in dz because of doubled depth bug
%                         % 1st index (x) = row index = y axis -- transposed -> x axis
%                         xlabel('x [m]')
%                         ylabel('z [m]')
%                 case 3
%                     imagesc(kgrid.x_vec, t_array*c0, squeeze(reflection_image(:,75,1:Nt))') % omit factor 1/2 in dz because of doubled depth bug
%                         xlabel('x [m]')
%                         ylabel('z [m]')
%             end
%                 title(['reconstructed image with c0 = ' num2str(c0) ', t0 correction = ' num2str(params.Nt_t0_correct)])
%                 cmap = colormap(gray);
%                 %cmap = flipud(cmap);    % flip colormap to make black = signal
%                 colormap(cmap);
%                 colorbar
%                 axis image
%                 drawnow


            %% find image peaks and FWHM for resolution measurements

%             kgrid.dt = params.dt;
% 
%             %-----
%             temp = load([dir_figures 'ROIstack_t0_' num2str(samples_t0_correct) '.mat']);
%             ROIstack = temp.ROIstack;
%             save([dir_figures 'ROIstack.mat'], 'ROIstack');
%             %-----
%             
%             threshold = 70;
%             [peaksInfo, ROI] = imagePeakFinder(reflection_image, c0, threshold, dir_figures, scanID, 'methodFWHM', '1DgaussianFitLat', 'reduceSensitivity', false);
% 
%             % concatenate peaksInfo array for all line scans
%             if(~exist('peaksInfoAll','var'))
%                 peaksInfoAll = peaksInfo;
%             else
%                 peaksInfoAll = cat(2,peaksInfoAll,peaksInfo);
%             end
% 
%             % stack ROI masks for all line scans
%             if(~exist([dir_figures 'ROIstack_t0_' num2str(samples_t0_correct) '.mat'],'file'))
%                 if(~exist('ROIstack','var'))
%                     ROIstack = cell(length(scanIDs));
%                 end
%                 ROIstack{scanID} = ROI;
%             end


            %% save figures

            % savefig(fig_data,[dir_figures 'autoplots\scan' num2str(scanID) '_sensor_data'], 'compact')
            % saveas(fig_data, [dir_figures 'autoplots\scan' num2str(scanID) '_sensor_data.jpg'])
            % 
            % savefig(fig_profile,[dir_figures 'autoplots\scan' num2str(scanID) '_profile'], 'compact')
            % saveas(fig_profile, [dir_figures 'autoplots\scan' num2str(scanID) '_profile.jpg'])
            % 
            % savefig(fig_image,[dir_figures 'autoplots\scan' num2str(scanID) '_image_marked'], 'compact')
            % saveas(fig_image, [dir_figures 'autoplots\scan' num2str(scanID) '_image_marked.jpg'])
            
%             savefig(fig_image,[dir_figures 'autoplots\scan' num2str(scanID) '_image_unmarked'], 'compact')
%             saveas(fig_image, [dir_figures 'autoplots\scan' num2str(scanID) '_image_unmarked.eps'])
            
            %% add to compound
            
            compound_image(1:size(reflection_image,1),1:size(reflection_image,2)) = compound_image(1:size(reflection_image,1),1:size(reflection_image,2)) + reflection_image;


            %% run multiple reconstructions in loop (end)

            % pause

        end     % of scanID / file_name loop

%         save( [dir_figures 'peaksInfoAll_x' file_name_xendings{idx_x} '_c' sprintf('%0.1f',c0) '_t0_' num2str(samples_t0_correct) '.mat'] , 'peaksInfoAll')
%         resoLat_contour_plot(peaksInfoAll, c0, samples_t0_correct, dir_figures)
%         clear peaksInfoAll
        
%         if(~exist([dir_figures 'ROIstack_t0_' num2str(samples_t0_correct) '.mat'],'file'))
%             save( [dir_figures 'ROIstack_t0_' num2str(samples_t0_correct) '.mat'] , 'ROIstack')
%         end
    
    end     % of c0 loop
    end     % of t0 loop
    
    
    %% plot compound_image
    
    global t_array Nt
    
    fig_image = figure;
    set(gcf,'Position',[100,100,800,450])
    imagesc(kgrid.x_vec, t_array*c0, compound_image(:,1:Nt)')
    xlabel('x [m]')
    ylabel('z [m]')
    title(['reconstructed image with c0 = ' num2str(c0) ', t0 correction = ' num2str(params.Nt_t0_correct)])
    cmap = colormap(gray);
    colormap(cmap);
    colorbar
    axis image
%     xlim([-0.0030,-0.0022])
%     ylim([0.0032,0.0040])
    
    
    %% USimagingAnalysis - movie of resoLat/countour plots at diff sound speeds
    
%     vidObj1 = VideoWriter([dir_figures '\vid_resoLat_c0.avi']);
%     vidObj2 = VideoWriter([dir_figures '\vid_resoLat_blur15_contour_c0.avi']);
%     
%     vidObj1.FrameRate = 5;
%     vidObj2.FrameRate = 5;
%     
%     open(vidObj1);
%     open(vidObj2);
% 
%     c0_sample = 1460:1:1500;
%     % c0_sample = 1482:0.2:1490;
%     % c0_sample = 1484:0.1:1487;
% 
%     for idx_c = 1:length(c0_sample)
% 
%         c0 = c0_sample(idx_c);
% 
%         figure(2*idx_c-1)
%         currFrame = getframe(gcf);
%         writeVideo(vidObj1,currFrame);
%         
%         figure(2*idx_c)
%         currFrame = getframe(gcf);
%         writeVideo(vidObj2,currFrame);
% 
%         pause(0.05)
% 
%     end
% 
%     close(vidObj1);
%     close(vidObj2);
%     
%     close all
    
end % of loop through 5 lines in x


%% save compound_image in .mat

filepath = 'D:\PROJECT\figures\_Matlab figs\USimaging\191126 resolution27umPlanar BK31[CNT] trolley scrambled fibre centralised parallel phantom\compound_image_x-05_c1489_t0-14.mat';

save(filepath, 'compound_image', 'kgrid', 't_array', 'params', 'c0', 'Nt')


%% load compound_image from .mat and imwrite to .tif

load(filepath)




%% t0 and c0 variation plots - grid overview

% t0_sample = -18:1:-12;
% c0_sample = 1484:1:1494;
% file_d='D:\PROJECT\figures\_Matlab figs\USimaging\191126 resolution27umPlanar BK31[CNT] trolley scrambled fibre centralised parallel phantom\x-05 c0_var t0_var gauss_offset\';
% 
% figure
% set(gcf,'Position',[70,0,1850,1000])
% 
% idx_subpl = 0;
% for idx_t = 1:length(t0_sample)
%     t0=t0_sample(idx_t);
%     for idx_c = 1:length(c0_sample)
%         c0 = c0_sample(idx_c);
%         
%         file_n = ['peaksInfoAll_x-05_c' num2str(c0) '.0_t0_' num2str(t0) '.mat'];
%         data = load([file_d file_n]);
%         peaksInfoAll = data.peaksInfoAll;
%         
%         peaksPosX    = peaksInfoAll(2,:)*1e3; % in mm
%         peaksPosZ    = peaksInfoAll(3,:)*1e3; % in mm
%         peaksResoLat = peaksInfoAll(4,:)*1e6; % in um
%         
%         [gridX,gridZ] = meshgrid(-11:0.1:11, 0:0.1:12.5);
%         resoLat = griddata(peaksPosX,peaksPosZ,peaksResoLat,gridX,gridZ,'linear');
%         
%         xreal = -11:0.1:11;
%         zreal = 0:0.1:12.5;
%         xbounds = 9:194;
%         zbounds = 9:124;
%         xreal = xreal(xbounds);
%         zreal = zreal(zbounds);
%         resoLat = resoLat(zbounds,xbounds);
%         
%         resoLatBlur = imgaussfilt(resoLat,15);
%         
%         idx_subpl = idx_subpl + 1;
%         subplot(7,11,idx_subpl)
%         contoursLatMajor = 45:5:140;
%         [Cmajor, hmajor]= contour(xreal, zreal, resoLatBlur, contoursLatMajor, 'LineWidth', 0.5);
%         clabel(Cmajor,hmajor, 'labelspacing', 700);
%         % clabel(C,'manual')
%         colormap(gray)
%         xlim([-8.3,8.3])
%         ylim([0.8,12.3])
%         
%     end
% end


%% sensor data image zoom in on parasitic signal to show no reverberations
% 
% file_dir = '..\data\imagingUS\191126\';
% file_name = 'resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_x-05_z0@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_29s48m13h_26-11-19_avg1_1D_raw.SGL';
% 
% [sensor_data, params] = loadSGL([file_dir file_name]);
% sensor_data = - sensor_data; % flip for trolley
% 
% time = linspace(1,params.Nt,params.Nt)*params.dt;
% time = time - time(14); % set t = 0 to 17 dt
% time = time * 1e6; % us
% 
% grid = kWaveGrid(params.Nx,params.dx,params.Ny,params.dy);
% 
% figure
% imagesc(grid.y_vec*1e3, time, sensor_data')
%     colormap(gray)
%     colorbar
%     xlabel('y axis [mm]')
%     ylabel('Time [\mus]')
%     xlim([2,7])
%     ylim([time(1),2.2])
%     set(gca,'FontName','Times New Roman')
%     set(gca,'FontSize',14)




