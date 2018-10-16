function viewSGL(dataSGL,params,t_0,slice_x,slice_y)
    
    Nx = params.Nx;     % size(dataSGL,1);
    Ny = params.Ny;     % size(dataSGL,2);
    Nt = size(dataSGL,3); % params.Nt;
    
    dx = params.dx;
    dy = params.dy;
    dt = params.dt;
    
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    time = (t_0+linspace(dt,Nt*dt,Nt))*1e6; % in us
    
    figure('pos',[250 150 800 500])
    
    %look at one slice through x = slice_x: y(t)
    subplot(2,2,1)
    imagesc(time,kgrid.x_vec*1e3,squeeze(dataSGL(slice_x,:,:)))
        title(['slice through x = ',num2str(kgrid.x_vec(slice_x)*1e3),' mm'])
        xlabel('time / \mu s')
        ylabel('y position / mm')
        colormap(gray)
        colorbar
        
    %look at one slice through y = slice_y: x(t)
    subplot(2,2,2)
    imagesc(time,kgrid.y_vec*1e3,squeeze(dataSGL(:,slice_y,:)))
        title(['slice through y = ',num2str(kgrid.y_vec(slice_y)*1e3),' mm'])
        xlabel('time / \mu s')
        ylabel('x position / mm')
        colormap(gray)
        colorbar
        
    %look at one time series at x/y = slice_x/y: signal-V(t)
    subplot(2,2,3)
    %data_V = dataSGL / 0.033 * 81.2 *1e-3 ;        % V-pressure conversion
    plot(time,squeeze(dataSGL(slice_x,slice_y,:)))
    %plot(squeeze(dataSGL(slice_x,slice_y,:)))
    %plot(data)
        title(['time series at x = ',num2str(kgrid.x_vec(slice_x)*1e3),' mm, y = ',num2str(kgrid.y_vec(slice_y)*1e3),' mm'])
        xlabel('time / \mu s')
        ylabel('signal amplitude / V')
        %ylim([-0.2,0.6])
    
    %look at MIP through time
    subplot(2,2,4)
    %dataSGL = -dataSGL;         % invert time series for positive peak
    peak_xy = max(dataSGL, [], 3);
    %dataSGL = dataSGL(50:110,30:110,:);
	%peak_xy = max(abs(dataPressure), [], 3);
    peak_max = max(max(peak_xy));
    peak_avg = mean(mean(peak_xy));
    disp(['max of peaks across sensor = ',num2str(peak_max),' V'])
    disp(['average of peaks across sensor = ',num2str(peak_avg),' V'])
    %[X,Y] = meshgrid(kgrid.x_vec*1e3,kgrid.y_vec*1e3);
    %surf(X',Y',peak_xy)
    imagesc(kgrid.x_vec*1e3,kgrid.y_vec*1e3,peak_xy)
    %imagesc(peak_xy)
        title('MIP through time')
        xlabel('y position / mm')
        ylabel('x position / mm')
        view([0,0,1])
        axis image
        caxis([0,peak_max])
        colorbar
    
end