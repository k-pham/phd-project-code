function [peaksInfo, ROI] = imagePeakFinder(reflection_image, c0, threshold, dir_figures, scanID, varargin)
% find peak positions and amplitude above given threshold, and FWHM at peak
% used for resolution measurements

    global kgrid t_array dt
    
    clear peaksSort peaksInfo
    
	%% set usage defaults
    num_req_input_variables = 3;
    methFWHM = 'peak';
    reduceSens = false;
    
    % replace with user defined values if provided
    if nargin < num_req_input_variables
        error('Incorrect number of inputs.');
    elseif ~isempty(varargin)
        for input_index = 1:2:length(varargin)
            switch varargin{input_index}
                case 'methodFWHM'
                    methFWHM = varargin{input_index + 1};
                case 'reduceSensitivity'
                    reduceSens = varargin{input_index + 1};
                otherwise
                    error('Unknown optional input.');
            end
        end
    end

    %% make spatial axes in m
    xaxis = kgrid.x_vec;
    zaxis = t_array*c0/2;

    %% determine image size in pixels to limit boundaries of search
    sizeX = size(reflection_image,1);
    sizeZ = size(reflection_image,2);

    %% determine necessary half-width in pixels to get HWHM of 200 um
    halfwidth = 200e-6;
    hwX = round(halfwidth/kgrid.dx);
    hwZ = round(halfwidth/(kgrid.dt*c0));
    
    %% threshold image and find clusters
    imageThresholdMask = reflection_image > threshold;
    
    if exist([dir_figures 'ROIstack.mat'],'file')
        data = load([dir_figures 'ROIstack.mat']);
        ROIstack = data.ROIstack;
        ROI = ROIstack{scanID};
    else
        ROI = roipoly;
    end
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
    peakRmRadius = 0.7e-3;    % in m
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
        
            target_rangeX = peakPosX-hwX : peakPosX+hwX;
            target_rangeZ = peakPosZ-hwZ : peakPosZ+hwZ;
            
            target_xaxis = xaxis(target_rangeX);
            target_zaxis = zaxis(target_rangeZ)*2; % factor of 2 from bug
            
            if ~reduceSens
                target_xprofile = reflection_image(target_rangeX,peakPosZ);
                target_zprofile = reflection_image(peakPosX,target_rangeZ);
            else
                target_xprofile = mean(reflection_image(target_rangeX,peakPosZ-1:peakPosZ+1),2);
                target_zprofile = mean(reflection_image(peakPosX-1:peakPosX+1,target_rangeZ),1);
            end
            target_xydata   = reflection_image(target_rangeX,target_rangeZ)';
            
            init_ampl = peakAmpl;
            init_xmu  = xaxis(peakPosX);
            init_lat  = 60e-6/(2*sqrt(log(2)));
            init_zmu  = zaxis(peakPosZ)*2;
            init_axi  = 40e-6/(2*sqrt(log(2)));
            
%             try
                switch methFWHM
                    % get resolution using fwhm around peak
                    case 'peak'
                        peakFWHMlateral = fwhm(target_xprofile,kgrid.dx,0);    % lateral resolution
                        peakFWHMaxial   = fwhm(target_zprofile,dt*c0,0);       % axial resolution
                                                                                                        % omit factor 1/2 bc of depth bug
                    % get resolution using 1d gaussian fit around peak
                    case '1DgaussianFit'
                        init_x = [init_ampl, init_xmu, init_lat];
                        fit_x = fit(target_xaxis, target_xprofile, 'gauss1', 'Start', init_x);
                        peakFWHMlateral = fit_x.c1*2*sqrt(log(2));
                        %target_xaxisUPS = target_xaxis(1):1e-6:target_xaxis(end);
                        %target_xfitUPS  = fit_x.a1*exp(-((target_xaxisUPS-fit_x.b1)/fit_x.c1).^2);
                        %figure; plot(target_xaxis,target_xprofile,'b+', target_xaxisUPS,target_xfitUPS,'r--')
                        
                        init_z = [init_ampl, init_zmu, init_axi];
                        fit_z = fit(target_zaxis', target_zprofile', 'gauss1', 'Start', init_z);
                        peakFWHMaxial   = fit_z.c1*2*sqrt(log(2));
                        %target_zaxisUPS = target_zaxis(1):1e-6:target_zaxis(end);
                        %target_zfitUPS  = fit_z.a1*exp(-((target_zaxisUPS-fit_z.b1)/fit_z.c1).^2);
                        %figure; plot(target_zaxis,target_zprofile,'b+', target_zaxisUPS,target_zfitUPS,'r--')
                    
                    % get lateral resolution using 1d gaussian fit around peak and axial resolution using fwhm around peak
                    case '1DgaussianFitLat'
                        init_x = [init_ampl, init_xmu, init_lat];
                        fit_x = fit(target_xaxis, target_xprofile, 'gauss1', 'Start', init_x);
                        peakFWHMlateral = fit_x.c1*2*sqrt(log(2));
                        peakFWHMaxial   = fwhm(target_zprofile,dt*c0,0);        % omit factor 1/2 bc of depth bug                     
                                                         
                    % get resolution using 2d gaussian fit
                    case '2DgaussianFit'
                        [target_xmesh, target_zmesh] = meshgrid(target_xaxis,target_zaxis);
                        target_mesh(:,:,1) = target_xmesh;
                        target_mesh(:,:,2) = target_zmesh;

                        coeffs_init = [init_ampl, init_xmu, init_lat, init_zmu, init_axi];
                        [coeffs] = lsqcurvefit(@fun2DGaussian,coeffs_init,target_mesh,target_xydata);
                        peakFWHMlateral = coeffs(3)*2*sqrt(2*log(2));
                        peakFWHMaxial   = coeffs(5)*2*sqrt(2*log(2));
                        
                    otherwise
                        error('Unknown FWHM method input.');
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
%             catch
%                 warning('some peaks did not work')
%             end % of try
        
        end % of if peak not too close to edge

    end % loop through all peaks in peaksSort

end