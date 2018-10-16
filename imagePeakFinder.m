function imagePeakFinder(reflection_image, kgrid, t_array, c0, threshold)
% find peak positions and amplitude above given threshold, and FWHM at peak
% used for resolution measurements

    % threshold image and find clusters
    imageThresholdMask = reflection_image > threshold;
    imageThresholdMask(:,1:50) = 0;         % exclude high ampl noise near source
%     imageThresholdMask(:,500:end) = 0;      % exclude high ampl back of sensor reflection
%     imageThresholdMask(1:100,:) = 0;        % exclude high ampl reflections off frame
    
%     imageThresholdMask(:,1:450) = 0;        % exclude high ampl noise near source
%     imageThresholdMask(:,740:820) = 0;      % exclude high ampl back of sensor
%     imageThresholdMask(:,900:end) = 0;      % exclude high ampl rest noise
%     imageThresholdMask(1:120,:) = 0;        % exclude high ampl reflections off frame
    
% 	imageThresholdMask(:,1:950) = 0;        % exclude high ampl noise near source
%     imageThresholdMask(:,1850:end) = 0;      % exclude high ampl rest noise
%     imageThresholdMask(1:100,:) = 0;        % exclude high ampl reflections off frame

    
%     ROI = roipoly;
%     imageThresholdMask(ROI) = 0;
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
    peakRmRadius = 0.5;    % in mm
    while index <= size(peaksSort,2)
        xpos = xaxis(peaksSort(2,:));
        zpos = zaxis(peaksSort(3,:))'*2;    % *2 factor from bug
        dist = sqrt( ( xpos-xpos(index) ).^2 + ...
                     ( zpos-zpos(index) ).^2 );
        peaksToKeep = dist > peakRmRadius;
        peaksToKeep(index) = 1;
        peaksSort = peaksSort(:,peaksToKeep);
        index = index +1;
    end
    
    % output peak info and mark on current plot
    for peakIndex = 1 : size(peaksSort,2)
        
        peakAmpl = peaksSort(1,peakIndex);
        peakPosX = peaksSort(2,peakIndex);
        peakPosZ = peaksSort(3,peakIndex);
%         peakFWHM = fwhm(reflection_image(peakPosX-30:peakPosX+30,peakPosZ),kgrid.dx,0); % lateral resolution
        peakFWHM = fwhm(reflection_image(peakPosX,peakPosZ-30:peakPosZ+30),kgrid.dx,0); % axial resolution
        
        figure(gcf)
        hold on
        %plot(peakPosX,peakPosZ,'r.')
        plot(xaxis(peakPosX),zaxis(peakPosZ)*2,'r+')    % *2 factor from bug
        text(xaxis(peakPosX),zaxis(peakPosZ)*2,num2str(peakFWHM*1e6,3))

        if peakIndex == 1 disp('cluster ##: peakAmpl [peakPosX,peakPosZ] fwhm'), end %#ok<SEPEX>
        disp(['cluster #',num2str(peakIndex),': ',num2str(peakAmpl),...
            ' [',num2str(xaxis(peakPosX)),', ',num2str(zaxis(peakPosZ)*2),'], ',num2str(peakFWHM*1e6,3),' um'])

%         pause
    end

end