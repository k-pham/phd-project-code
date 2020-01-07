function resoLat_contour_plot(peaksInfoAll, c0, t0_correct, dir_figures)

% peaksAmpl    = peaksInfoAll(1,:);
peaksPosX    = peaksInfoAll(2,:)*1e3; % in mm
peaksPosZ    = peaksInfoAll(3,:)*1e3; % in mm
peaksResoLat = peaksInfoAll(4,:)*1e6; % in um
% peaksResoAxi = peaksInfoAll(5,:)*1e6; % in um

[gridX,gridZ] = meshgrid(-11:0.1:11, 0:0.1:12.5);

resoLat = griddata(peaksPosX,peaksPosZ,peaksResoLat,gridX,gridZ,'linear');
% resoAxi = griddata(peaksPosX,peaksPosZ,peaksResoAxi,gridX,gridZ,'linear');
% amplitude = griddata(peaksPosX,peaksPosZ,peaksAmpl,gridX,gridZ,'linear');

fig_resoLat = figure;
set(gcf,'Position',[700,400,800,450])
imagesc(resoLat)
colorbar

    set(gca,'FontSize',13)
    title(['c0 = ' sprintf('%0.1f',c0) ])
    caxis([0,200])
    xlabel('x axis / 0.1 mm')
    ylabel('depth z / 0.1 mm')

% figure
% set(gcf,'Position',[100,100,800,450])
% imagesc(resoAxi)

xreal = -11:0.1:11;
zreal = 0:0.1:12.5;

% xbounds = 7:185;
xbounds = 9:194;
zbounds = 9:124;

xreal = xreal(xbounds);
zreal = zreal(zbounds);
resoLat = resoLat(zbounds,xbounds);

resoLat = fillmissing(resoLat,'nearest');
% resoAxi = fillmissing(resoAxi,'nearest');

% resoLat(isnan(resoLat)) = 70;
% resoAxi(isnan(resoAxi)) = 40;
% amplitude(isnan(amplitude)) = 0;

resoLatBlur = imgaussfilt(resoLat,15);
% resoAxiBlur = imgaussfilt(resoAxi,15);

contoursLat = [40:2:60,60:5:120,120:10:200,200:20:500];
fig_contour = figure;
set(gcf,'Position',[100,100,800,450])
[C, h]= contour(xreal, zreal, resoLatBlur, contoursLat, 'LineWidth', 2);
    clabel(C,h, 'labelspacing', 700);
    %clabel(C,'manual')
    colormap(gray)
    colorbar
    set(gca,'FontSize',13)
    title(['c0 = ' sprintf('%0.1f',c0) ' t0 = ' num2str(t0_correct) ])
    caxis([0,400])
    xlabel('x axis / mm')
    ylabel('depth z / mm')
drawnow

savefig(fig_contour,[dir_figures 'resoLat_contour_blur15_c' num2str(10*c0) '_t0_' num2str(t0_correct) '.fig'], 'compact')
saveas(fig_contour, [dir_figures 'resoLat_contour_blur15_c' num2str(10*c0) '_t0_' num2str(t0_correct) '.jpg'])

savefig(fig_resoLat,[dir_figures 'resoLat_c' num2str(10*c0) '_t0_' num2str(t0_correct) ], 'compact')
saveas(fig_resoLat, [dir_figures 'resoLat_c' num2str(10*c0) '_t0_' num2str(t0_correct) '.jpg'])

end