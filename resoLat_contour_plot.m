function resoLat_contour_plot(peaksInfoAll, c0)

% peaksAmpl    = peaksInfoAll(1,:);
peaksPosX    = peaksInfoAll(2,:)*1e3; % in mm
peaksPosZ    = peaksInfoAll(3,:)*1e3; % in mm
peaksResoLat = peaksInfoAll(4,:)*1e6; % in um
% peaksResoAxi = peaksInfoAll(5,:)*1e6; % in um

[gridX,gridZ] = meshgrid(-11:0.1:11, 0:0.1:12.5);

resoLat = griddata(peaksPosX,peaksPosZ,peaksResoLat,gridX,gridZ,'linear');
% resoAxi = griddata(peaksPosX,peaksPosZ,peaksResoAxi,gridX,gridZ,'linear');
% amplitude = griddata(peaksPosX,peaksPosZ,peaksAmpl,gridX,gridZ,'linear');

figure
set(gcf,'Position',[700,400,800,450])
imagesc(resoLat)

% figure
% set(gcf,'Position',[100,100,800,450])
% imagesc(resoAxi)

xreal = -11:0.1:11;
zreal = 0:0.1:12.5;

% xbounds = 7:185;
xbounds = 9:184;
zbounds = 15:121;

xreal = xreal(xbounds);
zreal = zreal(zbounds);
resoLat = resoLat(zbounds,xbounds);

% resoLat = fillmissing(resoLat,'nearest');
% resoAxi = fillmissing(resoAxi,'nearest');

% resoLat(isnan(resoLat)) = 70;
% resoAxi(isnan(resoAxi)) = 40;
% amplitude(isnan(amplitude)) = 0;

resoLatBlur = imgaussfilt(resoLat,15);
% resoAxiBlur = imgaussfilt(resoAxi,15);

contoursLat = [50:2:60,60:5:120];
figure
set(gcf,'Position',[100,100,800,450])
[C, h]= contour(xreal, zreal, resoLatBlur, contoursLat, 'LineWidth', 2);
    clabel(C,h, 'labelspacing', 700);
%     clabel(C,'manual')
    colormap(gray)
    colorbar
    title(['sound speed c0 = ' num2str(c0)])
%     xlim([4,15])
%     ylim([1.1,3.2])


end