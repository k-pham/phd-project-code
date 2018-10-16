function compareSGLmovie(dataSGL_cold, dataSGL_heat, params, t_0)

    Nx = params.Nx;     % size(dataSGL_cold,1);
    Ny = params.Ny;     % size(dataSGL_cold,2);
    Nt = params.Nt;     % size(dataSGL_cold,3);
    
    dx = params.dx;
    dy = params.dy;
    dt = params.dt;

    kgrid = kWaveGrid(Nx, dx, Ny, dy);
    time = (t_0+linspace(dt,Nt*dt,Nt))*1e6;
    difference = dataSGL_heat - dataSGL_cold;

    figure
    for slice_x = 1:Nx
        for slice_y = 1:Ny
            plot(time,squeeze(dataSGL_cold(slice_x,slice_y,:)),'b')
            hold on
            plot(time,squeeze(dataSGL_heat(slice_x,slice_y,:)),'r')
            plot(time,squeeze(difference(slice_x,slice_y,:)),'k--')
                title(['time series at x = ',num2str(kgrid.x_vec(slice_x)*1e3),' mm, y = ',num2str(kgrid.y_vec(slice_y)*1e3),' mm'])
                xlabel('time / \mu s')
                ylabel('signal amplitude / V')
                %legend('pulser triggering','laser-pulseGen-pulser triggering','difference')
                axis([8.3,10.3,-0.8,0.8])
            hold off
            drawnow
            pause(0.01)
        end
    end

end