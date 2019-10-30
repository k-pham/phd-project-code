function peaksInfo = imagePeakFinder(reflection_image, c0, threshold)
% find peak positions and amplitude above given threshold, and FWHM at peak
% used for resolution measurements

    global kgrid t_array dt
    
    clear peaksSort peaksInfo
    
    sizeX = size(reflection_image,1);
    sizeZ = size(reflection_image,2);
    
    % determine necessary half-width in pixels to get FWHM of 100 um
%     hwX = 15;
%     hwZ = 50;
    halfwidth = 100e-6;
    hwX = round(halfwidth/kgrid.dx);
    hwZ = round(halfwidth/(kgrid.dt*c0));

    % threshold image and find clusters
    imageThresholdMask = reflection_image > threshold;
%     imageThresholdMask(:,1:50) = 0;         % exclude high ampl noise near source
%     imageThresholdMask(:,500:end) = 0;      % exclude high ampl back of sensor reflection
%     imageThresholdMask(1:100,:) = 0;        % exclude high ampl reflections off frame
    
%     imageThresholdMask(:,1:450) = 0;        % exclude high ampl noise near source
%     imageThresholdMask(:,740:820) = 0;      % exclude high ampl back of sensor
%     imageThresholdMask(:,900:end) = 0;      % exclude high ampl rest noise
%     imageThresholdMask(1:120,:) = 0;        % exclude high ampl reflections off frame
    
%     imageThresholdMask(:,1:950) = 0;        % exclude high ampl noise near source
%     imageThresholdMask(:,1850:end) = 0;     % exclude high ampl rest noise
%     imageThresholdMask(1:50,:) = 0;         % exclude high ampl reflections off frame

%     imageThresholdMask(2000:end,:) = 0;     % exclude frame
%     imageThresholdMask(:,1:3500) = 0;
%     imageThresholdMask(:,4500:end) = 0;
    
    ROI = roipoly;
    imageThresholdMask(ROI' == 0) = 0;
    imageThresholdClusters = bwconncomp(imageThresholdMask);

    % make spatial axes in mm
    xaxis = kgrid.x_vec*1e3;
    zaxis = t_array*c0/2*1e3;
    
    % for each cluster go through pixels of cluster and find peak (ampl, pos, fwhm)
    for clusterIndex = 1 : imageThresholdClusters.NumObjects
        pixelIndexList = imageThresholdClusters.PixelIdxList{clusterIndex};
        clusterAmpl = reflection_image(pixelIndexList);
        peakAmpl = max(clusterAmpl);
        peakPixelIndex = pixelIndexList(clusterAmpl == peakAmpl);
        [peakPosX,peakPosZ] = ind2sub(size(reflection_image),peakPixelIndex);
        %[peakPosX,peakPosZ] = find(reflection_image == peakAmpl);    % also works
        
        peaks(1,clusterIndex) = peakAmpl;
        peaks(2,clusterIndex) = peakPosX;
        peaks(3,clusterIndex) = peakPosZ;
    end
    
    % sort peaks in order of amplitude
    [~, peakAmplSortIdx] = sort(peaks(1,:),'descend');
    peaksSort(1,:) = peaks(1,peakAmplSortIdx);
    peaksSort(2,:) = peaks(2,peakAmplSortIdx);
    peaksSort(3,:) = peaks(3,peakAmplSortIdx);
    
    % remove nearby side peaks
    index = 1;
    peakRmRadius = 0.4;    % in mm
    while index <= size(peaksSort,2)
        xpos = xaxis(peaksSort(2,:));
        zpos = zaxis(peaksSort(3,:))'*2;
        dist = sqrt( ( xpos-xpos(index) ).^2 + ...
                     ( zpos-zpos(index) ).^2 );
        peaksToKeep = dist > peakRmRadius;
        peaksToKeep(index) = 1;
        peaksSort = peaksSort(:,peaksToKeep);
        index = index +1;
    end
    
    % make peaksInfo array for output incl Ampl, XPos, ZPos in real space
    peaksInfo(1,:) = peaksSort(1,:);
    peaksInfo(2,:) = xaxis(peaksSort(2,:));
    peaksInfo(3,:) = zaxis(peaksSort(3,:))*2;
    
    % output peak info and mark on current plot
    for peakIndex = 1 : size(peaksSort,2)
        
        peakAmpl = peaksSort(1,peakIndex);
        peakPosX = peaksSort(2,peakIndex);
        peakPosZ = peaksSort(3,peakIndex);
        
        if (peakPosX - hwX > 0) && (peakPosX + hwX <= sizeX) && (peakPosZ - hwZ > 0) && (peakPosZ + hwZ <= sizeZ)
        
            try
                peakFWHMlateral = fwhm(reflection_image(peakPosX-hwX:peakPosX+hwX,peakPosZ),kgrid.dx,0);    % lateral resolution
                peakFWHMaxial   = fwhm(reflection_image(peakPosX,peakPosZ-hwZ:peakPosZ+hwZ),dt*c0,0);       % axial resolution
                                                                                                      % omit factor 1/2 bc of depth bug
                peaksInfo(4,peakIndex) = peakFWHMlateral;
                peaksInfo(5,peakIndex) = peakFWHMaxial;

                figure(gcf)
                hold on
                %plot(peakPosX,peakPosZ,'r.')
                plot(xaxis(peakPosX),zaxis(peakPosZ)*2,'r+')    % *2 factor from bug
                text(xaxis(peakPosX),zaxis(peakPosZ)*2,num2str(peakFWHMlateral*1e6,3),'HorizontalAlignment','left'  ,'VerticalAlignment','middle','Color','r')
                text(xaxis(peakPosX),zaxis(peakPosZ)*2,num2str(peakFWHMaxial*1e6,3)  ,'HorizontalAlignment','center','VerticalAlignment','top'   ,'Color','b')
                drawnow

                if peakIndex == 1 disp('cluster ##: peakAmpl [peakPosX,peakPosZ] fwhm: lateral / axial '), end %#ok<SEPEX>
                disp(['cluster #',num2str(peakIndex),': ',num2str(peakAmpl),...
                    ' [',num2str(xaxis(peakPosX)),', ',num2str(zaxis(peakPosZ)*2),'], ',...
                    num2str(peakFWHMlateral*1e6,3),' um, ',num2str(peakFWHMaxial*1e6,3),' um'])
            catch
                warning('some peaks did not work')
            end
        
        end
%         pause
    end

end