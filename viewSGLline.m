function viewSGLline(dataSGL,params,t_0,slice)

    Nx = params.Nx;     % size(dataSGL,1);
    Ny = params.Ny;     % size(dataSGL,1);
    Nt = params.Nt;     % size(dataSGL,2);
    
    dx = params.dx;
    dy = params.dy;
    dt = params.dt;
    
    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    timeAxis = (t_0+linspace(dt,Nt*dt,Nt))*1e6; % in us
    % decide whether line scan in y or in x to assign spatial axis etc
    if Nx == 1
        spaceAxis = kgrid.y_vec*1e3;
        positionString = ['y = ',num2str(kgrid.y_vec(slice)*1e3),' mm'];
    elseif Ny == 1
        spaceAxis = kgrid.x_vec*1e3;
        positionString = ['x = ',num2str(kgrid.x_vec(slice)*1e3),' mm'];
    end
    
    figure('pos',[80 50 600 800])


    % plot whole slice
    subplot(2,1,1)
    %data_V = dataSGL / 0.033 * 81.2 *1e-3 ;        % V-pressure conversion
    imagesc(timeAxis,spaceAxis,dataSGL(:,:))
        title('slice from line scan')
        xlabel('time / \mu s')
        ylabel('signal amplitude / V')
        colormap(gray)
        colorbar
        
    % plot single time series
    subplot(2,1,2)
    plot(timeAxis,dataSGL(slice,:))
        title(['time series at ',positionString])
        xlabel('time / \mu s')
        ylabel('signal amplitude / V')
%         ylim([-0.02,0.1])

end