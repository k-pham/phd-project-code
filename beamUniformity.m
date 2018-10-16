% spin off from coatingCharac.m specific to plotting beam profiles to
% assess uniformity

%% beam uniformity 180115
%assessed using moving excitation/diffuser-tube and single point
%interrogation and india ink as absorber, on all hard dielectric sensor
% 
% file_dir = '../data/coatingCharac/180115/';
% 
% position = 25:1:50;     % position of motor in mm w.r.t. homed position
% beamProfile = zeros(size(position));
% 
% figure
% for index = position
%     file_name = [ 'beamUniformity_' num2str(index) '.txt' ];
%     
%     fileID = fopen([file_dir file_name],'r');
%     formatSpec = '%f %f';
%     sizeData = [2 Inf];
%     data_timeseries = fscanf(fileID,formatSpec,sizeData);
%     fclose(fileID);
%     data_timeseries = data_timeseries';
% 
%     time = data_timeseries(:,1);    % in ns
%     vAC  = data_timeseries(:,2);    % in V
%     
%     realtime = 3 + time*1e-3;       % convert to real time in us (3=guess)
%     vAC  = -vAC;                    % invert for positive peak
%     
%     plot(realtime,vAC)
%     hold on
%         title('overlaid time series from beam profiling')
%         xlabel('time / us')
%         ylabel('AC signal / V')
%     
%     beamProfile(index-min(position)+1) = max(vAC);
%     
% end
% 
% figure
% plot(position, beamProfile)
%     title('beam profile of excitation through diffuser tube 1b')
%     xlabel('position of excitation+diffuser tube / mm')
%     ylabel('peak signal from Single Point Interrogation / V')
%     
%% beam uniformity & acquisition noise 180122
%assessed using moving excitation/diffuser-tube and single point
%interrogation and india ink as absorber, on all hard dielectric sensor
% 
% file_dir = '../data/coatingCharac/180122/';
% t_0 = 3;            % time delay after trigger in us
% 
% xpos = 34:1:50;     % position of motor in mm w.r.t. homed position
% ypos = 0:1:11;      % position of manual translation stage (inverted labelling direction to actual y axis)
% beamProfile = zeros(length(xpos),length(ypos));
% 
% % extract peaks for all time series and plot overlaid time series for each y
% for yindex = ypos
% %     figure
%     for xindex = xpos
%         file_name = [ 'beamUniformity_diffuser1b05_x' num2str(xindex) '_y' num2str(yindex) '.txt' ];
% 
%         fileID = fopen([file_dir file_name],'r');
%         formatSpec = '%f %f';
%         sizeData = [2 Inf];
%         data_timeseries = fscanf(fileID,formatSpec,sizeData);
%         fclose(fileID);
%         data_timeseries = data_timeseries';
% 
%         time = data_timeseries(:,1);    % in ns
%         vAC  = data_timeseries(:,2);    % in V
% 
%         realtime = t_0 + time*1e-3;     % convert to real time in us
%         vAC  = -vAC;                    % invert for positive peak
% 
% %         plot(realtime,vAC)
% %         hold on
% %             title(['overlaid time series from beam profiling with diffuser1b05 at y=' num2str(yindex)])
% %             xlabel('time / us')
% %             ylabel('AC signal / V')
% 
%         beamProfile(xindex-min(xpos)+1,yindex-min(ypos)+1) = max(vAC);
%     end
% end
% 
% %PLOT 1D beam profile for given y indeces - and make video
% figHand = figure;
% videoClass = VideoWriter('beamProfile_diffuser1b05.mp4','MPEG-4');
% videoClass.set('FrameRate',1)
% open(videoClass)
% for yindex = ypos
% plot(xpos, beamProfile(:,yindex-min(ypos)+1)')
% hold on
%     title(['1D beam profile of excitation through diffuser1b05 at y=' num2str(yindex)])
%     xlabel('x position of excitation+diffuser tube / mm')
%     ylabel('peak signal from Single Point Interrogation / V')
%     axis([34,50,0,0.14])
%     drawnow
%     frames = getframe(figHand);
%     writeVideo(videoClass,frames)
% end
% close(videoClass)
% 
% %PLOT 2D beam profile
% figure
% %[xmesh, ymesh] = meshgrid(xpos,ypos);
% surf(ypos, xpos, beamProfile)
% hold on
%     title('2D beam profile of excitation through diffuser1b05')
%     xlabel('y position of excitation+diffuser tube / mm')
%     ylabel('x position of excitation+diffuser tube / mm')
%     axis equal
%     %axis([0,11,34,50])
%     view([0,0,90])
%     colorbar
% 
% %ACQUISITION NOISE from repeated measurements of time series
% peakData = zeros(1,20);
% figure
% for xindex = 1:1:20
%     file_name = [ 'beamUniformity_diffuser1b05_x47_y5_noise' num2str(xindex) '.txt'];
% 
%     fileID = fopen([file_dir file_name],'r');
%     formatSpec = '%f %f';
%     sizeData = [2 Inf];
%     data_timeseries = fscanf(fileID,formatSpec,sizeData);
%     fclose(fileID);
%     data_timeseries = data_timeseries';
% 
%     time = data_timeseries(:,1);    % in ns
%     vAC  = data_timeseries(:,2);    % in V
%     
%     realtime = t_0 + time*1e-3;     % convert to real time in us
%     vAC  = -vAC;                    % invert for positive peak
% 
%     peakData(xindex) = max(vAC);
%     
%     plot(realtime,vAC)
%     hold on
%         title(['overlaid time series from beam profiling with diffuser1b05 at x=47, y=5 with peak = ' num2str(mean(peakData)) ' +/- ' num2str(std(peakData)) ])
%         xlabel('time / us')
%         ylabel('AC signal / V')
% 
% end

%% beam uniformity 180130
%assessed using moving excitation/diffuser-tube and single point
%interrogation and india ink as absorber, on all hard dielectric sensor
%beam profile acquired along cross

file_dir = '../data/coatingCharac/180130/';

xpos = 34:1:50;     % position of motor in mm w.r.t. homed position
ypos = 0:1:12;      % position of motor in mm w.r.t. homed position
% beamProfile_x       = zeros(1,length(xpos));
% beamProfile_y       = zeros(1,length(ypos));
beamProfile_x_avg10 = zeros(1,length(xpos));
beamProfile_y_avg10 = zeros(1,length(ypos));

for xindex = xpos
    yindex = 6;
    file_name = [ 'beamUniformity_diffuser1b05_x' num2str(xindex) '_y' num2str(yindex) '_avg10.txt' ];
    
    fileID = fopen([file_dir file_name],'r');
    formatSpec = '%f %f';
    sizeData = [2 Inf];
    data_timeseries = fscanf(fileID,formatSpec,sizeData);
    fclose(fileID);
    data_timeseries = data_timeseries';

    time = data_timeseries(:,1);    % in ns
    vAC  = data_timeseries(:,2);    % in V
    
    realtime = 5 + time*1e-3;       % convert to real time in us
    vAC  = -vAC;                    % invert for positive peak
    
    beamProfile_x_avg10(xindex-min(xpos)+1) = max(vAC);
    
end

for yindex = ypos
    xindex = 42;
    file_name = [ 'beamUniformity_diffuser1b05_x' num2str(xindex) '_y' num2str(yindex) '_avg10.txt' ];
    
    fileID = fopen([file_dir file_name],'r');
    formatSpec = '%f %f';
    sizeData = [2 Inf];
    data_timeseries = fscanf(fileID,formatSpec,sizeData);
    fclose(fileID);
    data_timeseries = data_timeseries';

    time = data_timeseries(:,1);    % in ns
    vAC  = data_timeseries(:,2);    % in V
    
    realtime = 5 + time*1e-3;       % convert to real time in us
    vAC  = -vAC;                    % invert for positive peak
    
    beamProfile_y_avg10(yindex-min(ypos)+1) = max(vAC);
    
end

figure
hold on
plot(xpos-35.5, beamProfile_x_avg10, 'b')
plot(ypos, beamProfile_y_avg10, 'r')
    title('beam profile of excitation through diffuser tube 1b')
    xlabel('position of excitation+diffuser tube along x or y axis / mm')
    ylabel('peak signal from Single Point Interrogation / V')
    legend('x axis','y axis')
    xlim([0,12])
    ylim([0,0.25])
    
