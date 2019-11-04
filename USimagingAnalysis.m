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

data = load('D:\PROJECT\figures\_Matlab figs\USimaging\190927 resolution27umPlanar BK31[CNT] trolley straight fibre/peaksInfoAll.mat');
peaksInfoAll = data.peaksInfoAll;

peaksAmpl    = peaksInfoAll(1,:);
peaksPosX    = peaksInfoAll(2,:);
peaksPosZ    = peaksInfoAll(3,:);
peaksResoLat = peaksInfoAll(4,:)*1e6; % in um
peaksResoAxi = peaksInfoAll(5,:)*1e6; % in um

[gridX,gridZ] = meshgrid(-11:0.1:11, 0:0.1:12.5);

resoLat = griddata(peaksPosX,peaksPosZ,peaksResoLat,gridX,gridZ,'linear');
resoAxi = griddata(peaksPosX,peaksPosZ,peaksResoAxi,gridX,gridZ,'linear');
amplitude = griddata(peaksPosX,peaksPosZ,peaksAmpl,gridX,gridZ,'linear');

resoLat = fillmissing(resoLat,'nearest');
resoAxi = fillmissing(resoAxi,'nearest');

resoLat(isnan(resoLat)) = 70;
resoAxi(isnan(resoAxi)) = 40;
amplitude(isnan(amplitude)) = 0;

resoLatBlur = imgaussfilt(resoLat,15);
resoAxiBlur = imgaussfilt(resoAxi,15);

figure
set(gcf,'Position',[100,100,800,450])
imagesc(resoLat)

figure
set(gcf,'Position',[100,100,800,450])
imagesc(resoAxi)

xreal = -11:0.1:11;
zreal = 0:0.1:12.5;

contoursLat = [50:2:60,60:5:120];
figure
set(gcf,'Position',[100,100,800,450])
[C, h]= contour(xreal, zreal, resoLatBlur, contoursLat, 'LineWidth', 2);
    clabel(C,h, 'labelspacing', 700);
%     clabel(C,'manual')
    colormap(gray)
    colorbar
%     xlim([4,15])
%     ylim([1.1,3.2])

contoursAxi = [30:5:45];figure
set(gcf,'Position',[100,100,800,450])
[C, h]= contour(xreal, zreal, resoAxiBlur, contoursAxi, 'LineWidth', 2);
    clabel(C,h, 'labelspacing', 700);
    colormap(gray)
    colorbar
    caxis([30,45])
    
