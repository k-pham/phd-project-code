% USimagingPhantoms.m:      parameters and file locations
% USimagingRecon.m:         run reconstruction, save .mat files, plot
% USimagingAnalysis.m:      image analysis


%% FWHM for LSF in resolution27
% 
% fwhm(reflection_image(:,201),kgrid.dx,1)
  
%% resolution27 - interpolate resolution and peak_amp from scattered data
% 
% x_vox = [ 402, 710, 1052, 1335, 1629, ...
%         181, 524, 848, 1163, 1443,    ...
%          308, 640, 945, 1250, 1552 ];
% z_vox = [ 200, 194, 211, 200, 201, ...
%         369, 365, 357, 365, 367,   ...
%          562, 539, 559, 550, 565 ];
% 
% resolution = [ 55.5, 56.6, 56.1, 55.9, 56.1, ...
%             74.7, 57.0, 57.5, 59.1, 59.1,    ...
%               73.5, 61.7, 64.1, 66.4, 68.5 ];
% 
% peak_amp = [ 110, 534, 955, 565, 42, ...
%            43, 179, 587, 682, 268,  ...
%             47, 397, 548, 498, 78 ];
% 
% [xq_vox,zq_vox] = meshgrid(0:1:1924, 0:1:750);
% 
% resolution_q = griddata(x_vox,z_vox,resolution,xq_vox,zq_vox,'cubic');
% peak_amp_q = griddata(x_vox,z_vox,peak_amp,xq_vox,zq_vox,'cubic');
% 
% resolution_q(isnan(resolution_q)) = 0;
% peak_amp_q(isnan(peak_amp_q)) = 0;

%% heat map of peak_amp for FOV
% 
% xreal = xq_vox * kgrid.dx *1e3;
% zreal = zq_vox * c0*dt *1e3;
% 
% figure(11); clf(11)
% imagesc(xreal(1,:),zreal(:,1),peak_amp_q)
% colormap(hot)
% colorbar
% xlim([4,15])
% ylim([1.17,3.2])

%% contour plot of resolution
% 
% contours = [ 56:1:60, 62:2:74 ];
% 
% xreal = xq_vox * kgrid.dx *1e3;
% zreal = zq_vox * c0*dt *1e3;
% 
% figure(13); clf(13)
% [C h]= contour(xreal, zreal, resolution_q, contours, 'LineWidth', 2);
% %clabel(C,h, 'labelspacing', 700);
% %clabel(C,'manual')
% colormap(gray)
% colorbar
% xlim([4,15])
% ylim([1.1,3.2])

%% layerGelWax - plot and fit exponential SNR
% 
% SNR = [824 348 291 223 220 145 ];
% z = [0.95 6.13 6.55 8.6 9.7 14.8];
% SNR = SNR/40;
% z = z+6;
% f = fit(z,SNR,'exp1');
% plot(f,z,SNR)


%% resolution analysis automated -- outsourced to imagePeakFinder.m

threshold = 100;
imagePeakFinder(reflection_image, kgrid, t_array, c0, threshold)

% imageThresholdMask = reflection_image > threshold;
% imageThresholdMask(:,1:180) = 0;    % exclude high ampl noise near source
% imageThresholdMask(:,260:end) = 0;    % exclude high ampl back of sensor reflection
% imageThresholdClusters = bwconncomp(imageThresholdMask);
% 
% xaxis = kgrid.x_vec*1e3;
% zaxis = t_array*c0/2*1e3;
% 
% % for each cluster go through pixels of cluster and find peak (ampl, pos)
% for clusterIndex = 1 : imageThresholdClusters.NumObjects
%     pixelIndexList = imageThresholdClusters.PixelIdxList{clusterIndex};
%     clusterAmpl = reflection_image(pixelIndexList);
%     peakAmpl = max(clusterAmpl);
%     peakPixelIndex = pixelIndexList(clusterAmpl == peakAmpl);
%     [peakPosX,peakPosZ] = ind2sub(size(reflection_image),peakPixelIndex);
%     %[peakPosX,peakPosZ] = find(reflection_image == peakAmpl);    % also works
%     
%     figure(gcf)
%     hold on
%     %plot(peakPosX,peakPosZ,'r.')
%     plot(xaxis(peakPosX),zaxis(peakPosZ)*2,'r.')    % *2 factor from bug
%     text(xaxis(peakPosX),zaxis(peakPosZ)*2,num2str(clusterIndex))
%     
%     peakFWHM = fwhm(reflection_image(peakPosX-50:peakPosX+50,peakPosZ),kgrid.dx,0);
%     
%     if clusterIndex == 1 disp('cluster ##: peakAmpl [peakPosX,peakPosZ] FWHM'), end %#ok<SEPEX>
%     disp(['cluster #',num2str(clusterIndex),': ',num2str(peakAmpl),...
%         ' [',num2str(xaxis(peakPosX)),', ',num2str(zaxis(peakPosZ)*2),'], ',num2str(peakFWHM*1e6),' um'])
%     
%     pause
% end


%% resolution27umPlanar contour plot
% use peaksInfo array

% data = load('D:\PROJECT\figures\_Matlab figs\USimaging\191029 resolution27umPlanar BK31[CNT] trolley scrambled fibre\peaksInfoAll.mat');
% data = load('D:\PROJECT\figures\_Matlab figs\USimaging\191031 resolution27umPlanar BK31[CNT] trolley scrambled fibre central phantom\peaksInfoAll.mat');
% data = load('D:\PROJECT\figures\_Matlab figs\USimaging\191031 resolution27umPlanar BK31[CNT] trolley scrambled fibre central phantom\meth = 1DgaussianFitLat\peaksInfoAll.mat');
% data = load('D:\PROJECT\figures\_Matlab figs\USimaging\191126 resolution27umPlanar BK31[CNT] trolley scrambled fibre centralised parallel phantom\x-1\peaksInfoAll_x-1_c1488.0.mat');
% data = load('D:\PROJECT\figures\_Matlab figs\USimaging\191126 resolution27umPlanar BK31[CNT] trolley scrambled fibre centralised parallel phantom\x-05\peaksInfoAll_x-05_c1488.0.mat');

file_dir = 'D:\PROJECT\figures\_Matlab figs\USimaging\191126 resolution27umPlanar BK31[CNT] trolley scrambled fibre centralised parallel phantom\x-05 c0_var t0_var\';
file_name = 'peaksInfoAll_x-05_c1488.0_t0_-13.mat';

data = load([file_dir file_name]);
peaksInfoAll = data.peaksInfoAll;

% peaksAmpl    = peaksInfoAll(1,:);
peaksPosX    = peaksInfoAll(2,:)*1e3; % in mm
peaksPosZ    = peaksInfoAll(3,:)*1e3; % in mm
peaksResoLat = peaksInfoAll(4,:)*1e6; % in um
% peaksResoAxi = peaksInfoAll(5,:)*1e6; % in um

[gridX,gridZ] = meshgrid(-11:0.1:11, 0:0.1:12.5);

resoLat = griddata(peaksPosX,peaksPosZ,peaksResoLat,gridX,gridZ,'linear');
% resoAxi = griddata(peaksPosX,peaksPosZ,peaksResoAxi,gridX,gridZ,'linear');
% amplitude = griddata(peaksPosX,peaksPosZ,peaksAmpl,gridX,gridZ,'linear');

% resoLat = fillmissing(resoLat,'nearest');
% resoAxi = fillmissing(resoAxi,'nearest');

% resoLat(isnan(resoLat)) = 70;
% resoAxi(isnan(resoAxi)) = 40;
% amplitude(isnan(amplitude)) = 0;

figure
set(gcf,'Position',[100,100,800,450])
imagesc(resoLat)

% figure
% set(gcf,'Position',[100,100,800,450])
% imagesc(resoAxi)

pause

xreal = -11:0.1:11;
zreal = 0:0.1:12.5;

% xbounds = 7:185;
% xbounds = 9:184;
xbounds = 9:194;
% zbounds = 15:121;
zbounds = 9:124;

xreal = xreal(xbounds);
zreal = zreal(zbounds);
resoLat = resoLat(zbounds,xbounds);

figure
set(gcf,'Position',[100,100,800,450])
imagesc(resoLat)

pause

resoLatBlur = imgaussfilt(resoLat,15);
% resoAxiBlur = imgaussfilt(resoAxi,15);

figure
set(gcf,'Position',[100,100,800,450])
    contoursLatMajor = 45:5:140;
    contoursLatMinor = setdiff(45:1:140,contoursLatMajor);
    hold on
    [Cmajor, hmajor]= contour(xreal, zreal, resoLatBlur, contoursLatMajor, 'LineWidth', 1);
    [Cminor, hminor] = contour(xreal,zreal,resoLatBlur, contoursLatMinor, 'LineWidth', 1, 'LineStyle', ':');
    hold off
    colormap(gray)
    xlim([-8.3,8.3])
    ylim([0.8,12.3])
    clabel(Cmajor,hmajor, 'labelspacing', 500);
%     clabel(Cmajor,'manual')
    
% contoursAxi = 30:5:45;
% figure
% set(gcf,'Position',[100,100,800,450])
% [Cmajor, hmajor]= contour(xreal, zreal, resoAxiBlur, contoursAxi, 'LineWidth', 2);
%     clabel(Cmajor,hmajor, 'labelspacing', 700);
%     colormap(gray)
%     colorbar
%     caxis([30,45])
    

%% movie of plots of depth profiles at diff latitudes

vidObj = VideoWriter('depthprofiles_latitude_blur25.avi');
open(vidObj);

figure
hold on
set(gca,'FontSize',13)
xlabel('depth z / mm')
ylabel('lateral resolution / \mum')
axis([0,12.5,50,180])
for latitude = 1:length(xreal)
    plot(zreal,resoLatBlur(:,latitude))
    title(['x position: '  sprintf('%0.1f',xreal(latitude)) ' mm'])
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    drawnow
    pause(0.05)
end

close(vidObj);


%% movie of resoLat/countour plots at diff sound speeds

vidObj = VideoWriter([dir_figures '\vid_resoLat_c0_fine2.avi']);
vidObj = VideoWriter([dir_figures '\vid_resoLat_blur15_contour_c0_fine2.avi']);
vidObj.FrameRate = 5;
open(vidObj);

% c0_sample = 1460:1:1500;
% c0_sample = 1482:0.2:1490;
c0_sample = 1484:0.1:1487;

for idx_c = 1:length(c0_sample)
    
    c0 = c0_sample(idx_c);
    
%     figure(2*idx_c-1)
    figure(2*idx_c)
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
    pause(0.05)
    
end

close(vidObj);






