
clearvars

dx = 5e-6;
sp = 1500;
Nx = 1512;
Ny = 384;
dy = dx;
kgrid = kWaveGrid(Nx, dx, Ny, dy);
clear source;
medium.sound_speed_ref =  sp;
medium.sound_speed = sp*ones(Nx,Ny);
medium.density = 1000*ones(Nx,Ny);
dt = kgrid.dt;
kgrid.t_array = 0:3e-10:10000*4e-10;

    ab = 15;

[x,y] = meshgrid(1:Ny,1:Nx);

i = 1;
for j = 1:7
       
    for i = 1:2
                      dist = sqrt((y - (100 + 200*(j-1) + 100*(i-1))).^2+(x-(100+(ab+6))  - 200*(i-1)).^2)*dx;
                       
                      medium.sound_speed(dist<13.5e-6) = 5200;
                      medium.density(dist<13.5e-6) = 6000;         
    end
end
   
        
    vec = 1:Nx';
    
    gauss = exp(-(vec-Nx/2).^2/(2*(0.6*(Nx/2)^2)));

    medium.sound_speed(:,ab:ab+6) = 1000;
    medium.density(:,ab:ab+6) = 1200;

    source.p0 = zeros(Nx,Ny);    
    % preallocate an empty pressure source matrix
    source.p0(:,15) = gauss;

sensor.mask = zeros(kgrid.Nx, kgrid.Ny);
sensor.mask(:,14) = 1;
sensor.mask(:,360) = 1;

sensor.record = {'p_max_all','p'};
 
Time_Reversed = kspaceFirstOrder2D(kgrid, medium, source, sensor,'PMLSize',[10,10],  'PlotScale', max(source.p0(:))*[-1, 1]);
P_Data = reshape(Time_Reversed.p,[Nx,2,size(Time_Reversed.p,2)]);

P_Data = permute(P_Data,[1,3,2]);
P_Data(:,1:200,:) = 0;

for j = 1:size(P_Data,1)
    P_Data2(j,:) = applyFilter(P_Data(j,:,1),1/3e-10,50e6,'LowPass','TransitionWidth',0.03);
end


[reflection_image, samples_total, dt, t_array, kgrid] = reconstruct3dUSimage(P_Data2(1:4:end,1:10:end,1), 4*dx, 0, 0, 0, 1460,3e-10*10);



ref = reflection_image(:,1:end/2);
figure; imagesc(-ref); axis image; colormap(gray);
