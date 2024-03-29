function p = kspacePlaneRecon_US(p, dy, dz, dt, c, varargin)


%KSPACEPLANERECON_US 3D planar FFT reconstruction for plane wave ultrasound
%                    reflection imaging 
%
% DESCRIPTION:
%     adapted from kspacePlaneRecon (which is for photoacoustic imaging).
%
% USAGE:
%     p_xyz = kspacePlaneRecon_US(p_tyz, dy, dz, dt, c)
%     p_xyz = kspacePlaneRecon_US(p_tyz, dy, dz, dt, c, ...)
%
% INPUTS:
%     p_tyz       - pressure time-series recorded over an evenly spaced
%                   array of sensor points on a line (indexed as t, y, z)
%     dy, dz      - spatial step [m]
%     dt          - time step [s]
%     c           - acoustically-homogeneous sound speed [m/s]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings.
%
%     'DataOrder' - String input which sets the data order (default =
%                   'tyz'). Valid inputs are 'tyz' and 'yzt'.
%     'Interp'    - String input controlling the interpolation method
%                   used by interp3 in the reconstruction (default =
%                   '*nearest').
%     'Plot'      - Boolean controlling whether a plot of the reconstructed
%                   estimate of the initial acoustic pressure distribution
%                   is produced (default = false). 
%     'PosCond'   - Boolean controlling whether a positivity condition is
%                   enforced on the reconstructed estimate of the initial
%                   acoustic pressure distribution (default = false).
%
% OUTPUTS:
%     p_xyz       - image (indexed as x, y, z) 
%
% ABOUT:
%     author      - Bradley Treeby and Ben Cox
%     date        - 2nd July 2009
%     last update - 21st June 2017
%       


% avoid interp3 error in Matlab R2012b
matlab_ver = ver('matlab'); 
if strcmp(matlab_ver.Release, '(R2012b)')
    error('This function will not run in Matlab R2012b because of a bug in Matlab''s function interp3. It does work in prior and subsequent versions of Matlab.')
end

% start timer
tic;

% define defaults
num_req_inputs = 5;
data_order = 'tyz';
interp_method = '*nearest';
plot_recon = false;
positivity_cond = false;

% replace with user defined values if provided
if nargin < num_req_inputs
    error('Incorrect number of inputs.');
elseif ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}          
            case 'DataOrder'
                data_order = varargin{input_index + 1};
                if ~strcmp(data_order, 'tyz') && ~strcmp(data_order, 'yzt')
                    error('Unknown setting for optional input DataOrder.');
                end
            case 'Interp'
                interp_method = varargin{input_index + 1};            
            case 'Plot'
                plot_recon = varargin{input_index + 1};
            case 'PlotRecon'
                plot_recon = varargin{input_index + 1};
            case 'PosCond'
                positivity_cond = varargin{input_index + 1};  
            otherwise
                error(['Unknown optional input ' varargin{input_index} '.']);
        end
    end
end

% reorder the data to p(t, y, z) if needed
if strcmp(data_order, 'yzt')
    p = permute(p, [3, 1, 2]);
end

% mirror the time domain data about t = 0 to allow the cosine transform in
% the t direction to be computed using an FFT
p = [flipdim(p, 1); p(2:end, :, :)];

% extract the size of mirrored input data
[Nt, Ny, Nz] = size(p);

% update command line status
disp('Running k-Wave planar reconstruction...');
disp(['  grid size: ' num2str((Nt + 1) / 2) ' by ' num2str(Ny) ' by ' num2str(Nz) ' grid points']);
disp(['  interpolation mode: ' interp_method]);

% create a computational grid that is evenly spaced in w, ky, and kz, where 
% Nx = Nt and dx = dt*c
kgrid = kWaveGrid(Nt, dt * c, Ny, dy, Nz, dz);

% from the grid for kx, create a computational grid for w using the
% relation dx = dt*c; this represents the initial sampling of p(w, ky, kz)
w = c .* kgrid.kx;

% remap the computational grid for kx onto w using the dispersion
% relation w/c = (kx^2 + ky^2 + kz^2)^1/2. This gives an w grid that is
% evenly spaced in kx. This is used for the interpolation from p(w, ky, kz)
% to p(kx, ky, kz). Only real w is taken to force kx (and thus x) to be
% symmetrical about 0 after the interpolation. 
%w_new = c .* kgrid.k;                      % photoacoustics
w_new = c .* kgrid.k.^2 ./ (2 * kgrid.kx) ; % reflection US
w_new(kgrid.kx==0) = 0;                     % reflection US

% calculate the scaling factor using the value of kx, where
% kx = sqrt( (w/c).^2 - kgrid.ky.^2 - kgrid.kz.^2 ) and then manually
% replacing the DC value (in the PA case) with its limit otherwise NaN results
% sf = c.^2 .* sqrt( (w ./ c).^2 - kgrid.ky.^2 - kgrid.kz.^2) ./ (2 .* w); %photoacoustic
% sf(w == 0 & kgrid.ky == 0 & kgrid.kz == 0) = c ./ 2;                     %photoacoustic
sf = kgrid.k + sqrt(kgrid.k.^2 - kgrid.ky.^2 - kgrid.kz.^2);               % planewave US

% compute the FFT of the input data p(t, y, z) to yield p(w, ky, kz) and
% scale
p = sf .* fftshift(fftn(ifftshift(p)));

% remove unused variables
clear sf;

% exclude the inhomogeneous part of the wave
p(abs(w) < (c * sqrt(kgrid.ky.^2 + kgrid.kz.^2))) = 0;

% compute the interpolation from p(w, ky, kz) to p(kx, ky, kz); for a
% matrix indexed as [M, N, P], the axis variables must be given in the
% order N, M, P
p = interp3(kgrid.ky, w, kgrid.kz, p, kgrid.ky, w_new, kgrid.kz, interp_method);

% remove unused variables
clear kgrid w;

% set values outside the interpolation range to zero
p(isnan(p)) = 0;

% compute the inverse FFT of p(kx, ky, kz) to yield p(x, y, z)
p = real(fftshift(ifftn(ifftshift(p))));

% remove the left part of the mirrored data which corresponds to the
% negative part of the mirrored time data
p = p( ((Nt + 1) / 2):Nt, :, :);

% correct the scaling - the forward FFT is computed with a spacing of dt
% and the reverse requires a spacing of dz = dt*c, the reconstruction
% assumes that p0 is symmetrical about z, and only half the plane collects
% data (first approximation to correcting the limited view problem) (p_zxy)
p = 2 * 2 * p ./ c;

% enfore positivity condition (p_zxy)
if positivity_cond
    disp('  applying positivity condition...');
    p(p < 0) = 0;
end

% update command line status
disp(['  computation completed in ' scaleTime(toc)]);

% plot the reconstruction
if plot_recon
    
    % allocate axis dimensions
    x_axis = [0, (Nt / 2) * dt * c];
    y_axis = [0, Ny * dy];
    z_axis = [0, Nz * dz];
    
    % select suitable axis scaling factor
    [~, scale, prefix] = scaleSI(max([x_axis(end), y_axis(end), z_axis(end)])); 
    
    % select suitable plot scaling factor
    plot_scale = max(p(:));
    
    % create the figures
    figure;
    subplot(2, 2, 1);
    imagesc(y_axis * scale, x_axis * scale, squeeze(p(:, :, round(end / 2))), [-plot_scale, plot_scale]);
    axis image;
    title('x-y plane');
    
    subplot(2, 2, 2);
    imagesc(z_axis * scale, x_axis * scale, squeeze(p(:, round(end / 2), :)), [-plot_scale, plot_scale]);
    axis image;
    xlabel(['(All axes in ' prefix 'm)']);
    title('x-z plane');
    
    subplot(2, 2, 3);
    imagesc(z_axis * scale, y_axis * scale, squeeze(p(round(end / 2), :, :)), [-plot_scale, plot_scale]);
    axis image;
    title('y-z plane');
    colormap(getColorMap);
    
end 