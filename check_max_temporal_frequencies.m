% script to check maximum temporal frequencies involved in
% kspaceLineRecon_US reconstruction.

% set-up acquisition and array parameters
Nt = 2500;
dt = 4e-9;
c = 1500;
Ny = 100;
dy = 150e-6;

disp(['For a ' num2str(Nt*dt*1e6) '-us acquisition, and ' num2str(Ny*dy*1e3) '-mm aperture we have:'])

% max detectable temporal frequency
kgrid_rec = kWaveGrid(Nt, dt*c, Ny, dy);
w = c .* kgrid_rec.kx;
disp(['max detectable frequency: ' num2str(max(w,[],'all')/1e6) ' MHz'])

% max temporal frequency that uniform kz grid is mapped to
% if factor 1/2 included in object kgrid
kgrid_obj = kWaveGrid(Nt, dt*c/2, Ny, dy);
w_new = c .* kgrid_obj.k.^2 ./ (2 * kgrid_obj.kx) ;
w_new(kgrid_obj.kx==0) = 0;
disp(['max mapping frequency with factor 1/2 in object kgrid: ' num2str(max(w_new,[],'all')/1e6) ' MHz'])

% if factor 1/2 omitted in object kgrid
kgrid_obj = kWaveGrid(Nt, dt*c, Ny, dy);
w_new = c .* kgrid_obj.k.^2 ./ (2 * kgrid_obj.kx) ;
w_new(kgrid_obj.kx==0) = 0;
disp(['max mapping frequency without factor 1/2 in object kgrid: ' num2str(max(w_new,[],'all')/1e6) ' MHz'])
