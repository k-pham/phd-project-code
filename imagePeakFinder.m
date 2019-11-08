function peaksInfo = imagePeakFinder(reflection_image, c0, threshold, varargin)
% find peak positions and amplitude above given threshold, and FWHM at peak
% used for resolution measurements

    global kgrid t_array dt
    
    clear peaksSort peaksInfo
    
	%% set usage defaults
    num_req_input_variables = 3;
    toFitGaussian = false;

    % replace with user defined values if provided
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'FitGaussian'
                    toFitGaussian = varargin{input_index + 1};
                otherwise
                    error('Unknown optional input.');
            end
        end
    end

    %% make spatial axes in mm
    xaxis = kgrid.x_vec*1e3;
    zaxis = t_array*c0/2*1e3;

    %% determine image size in pixels to limit boundaries of search
    sizeX = size(reflection_image,1);
    sizeZ = size(reflection_image,2);

    %% determine necessary half-width in pixels to get HWHM of 100 um
    halfwidth = 100e-6;
    hwX = round(halfwidth/kgrid.dx);
    hwZ = round(halfwidth/(kgrid.dt*c0));
    
    %% threshold image and find clusters
    imageThresholdMask = reflection_image > threshold;
    
    ROI = roipoly;
    imageThresholdMask(ROI' == 0) = 0;
    imageThresholdClusters = bwconncomp(imageThresholdMask);

    %% for each cluster go through pixels of cluster and find peak (ampl, pos, fwhm)
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
    
    %% sort peaks in order of amplitude
    [~, peakAmplSortIdx] = sort(peaks(1,:),'descend');
    peaksSort(1,:) = peaks(1,peakAmplSortIdx);
    peaksSort(2,:) = peaks(2,peakAmplSortIdx);
    peaksSort(3,:) = peaks(3,peakAmplSortIdx);
    
    %% remove nearby side peaks starting from high ampl peaks
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
    
    %% make peaksInfo array for output incl Ampl, XPos, ZPos in real space
    peaksInfo(1,:) = peaksSort(1,:);
    peaksInfo(2,:) = xaxis(peaksSort(2,:));
    peaksInfo(3,:) = zaxis(peaksSort(3,:))*2;
    
    %% get resolution, output peak info and mark on current plot
    for peakIndex = 1 : size(peaksSort,2)
        
        peakAmpl = peaksSort(1,peakIndex);
        peakPosX = peaksSort(2,peakIndex);
        peakPosZ = peaksSort(3,peakIndex);
        
        % check peak not too close to edge
        if (peakPosX - hwX > 0) && (peakPosX + hwX <= sizeX) && (peakPosZ - hwZ > 0) && (peakPosZ + hwZ <= sizeZ)
        
            try
                switch toFitGaussian
                    % get resolution using fwhm around peak
                    case false
                        peakFWHMlateral = fwhm(reflection_image(peakPosX-hwX:peakPosX+hwX,peakPosZ),kgrid.dx,0);    % lateral resolution
                        peakFWHMaxial   = fwhm(reflection_image(peakPosX,peakPosZ-hwZ:peakPosZ+hwZ),dt*c0,0);       % axial resolution
                                                                                                           % omit factor 1/2 bc of depth bug
                    % get resolution using 1d gaussian fit around peak
                    case '1D'
                        continue

                    % get resolution using 2d gaussian fit
                    case '2D'
                        target_rangeX = peakPosX-hwX : peakPosX+hwX;
                        target_rangeZ = peakPosZ-hwZ : peakPosZ+hwZ;

                        target_xaxis = xaxis(target_rangeX);
                        target_zaxis = zaxis(target_rangeZ)*2;
                        [target_xmesh, target_zmesh] = meshgrid(target_xaxis,target_zaxis);
                        target_mesh(:,:,1) = target_xmesh;
                        target_mesh(:,:,2) = target_zmesh;
                        target_data  = reflection_image(target_rangeX,target_rangeZ)';

                        init_ampl = peakAmpl;
                        init_x    = xaxis(peakPosX);
                        init_lat  = 70e-6/(2*sqrt(2*log(2)));
                        init_z    = zaxis(peakPosZ)*2;
                        init_axi  = 40e-6/(2*sqrt(2*log(2)));
                        coeffs_init = [init_ampl, init_x, init_lat, init_z, init_axi];
                        [coeffs] = lsqcurvefit(@fun2DGaussian,coeffs_init,target_mesh,target_data);
                        peakFWHMlateral = coeffs(3)*2*sqrt(2*log(2));
                        peakFWHMaxial   = coeffs(5)*2*sqrt(2*log(2));
                    otherwise
                        error('Unknown optional input.');
                end

                %% save resolution in peaksInfo and label plot in gcf              
                peaksInfo(4,peakIndex) = peakFWHMlateral;
                peaksInfo(5,peakIndex) = peakFWHMaxial;

                figure(gcf)
                hold on
                %plot(peakPosX,peakPosZ,'r.')
                plot(xaxis(peakPosX),zaxis(peakPosZ)*2,'r+')    % *2 factor from bug
                text(xaxis(peakPosX),zaxis(peakPosZ)*2,[' ',num2str(peakFWHMlateral*1e6,3)],'HorizontalAlignment','left'  ,'VerticalAlignment','middle','Color','r')
                text(xaxis(peakPosX),zaxis(peakPosZ)*2,num2str(peakFWHMaxial*1e6,3)  ,'HorizontalAlignment','center','VerticalAlignment','top'   ,'Color','b')
                drawnow

                if peakIndex == 1 disp('cluster ##: peakAmpl [peakPosX,peakPosZ] fwhm: lateral / axial '), end %#ok<SEPEX>
                disp(['cluster #',num2str(peakIndex),': ',num2str(peakAmpl),...
                    ' [',num2str(xaxis(peakPosX)),', ',num2str(zaxis(peakPosZ)*2),'], ',...
                    num2str(peakFWHMlateral*1e6,3),' um, ',num2str(peakFWHMaxial*1e6,3),' um'])
            catch
                warning('some peaks did not work')
            end % of try
        
        end % of if peak not too close to edge

    end % loop through all peaks in peaksSort

end