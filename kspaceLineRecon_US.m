function p = kspaceLineRecon_US(p, dy, dt, c, varargin)
%KSPACELINERECON_US 2D linear FFT reconstruction for plane wave ultrasound
%                   reflection imaging 
%
% DESCRIPTION
%     based on kspaceLineRecon
% 
% USAGE:
%     p_xy = kspaceLineRecon_US(p_ty, dy, dt, c)
%     p_xy = kspaceLineRecon_US(p_ty, dy, dt, c, ...)
%
% INPUTS:
%     p_ty        - pressure time-series recorded over an evenly spaced
%                   array of sensor points on a line (indexed as t, y)
%     dy          - spatial step [m]
%     dt          - time step [s]
%     c           - acoustically-homogeneous sound speed [m/s]
%
% OPTIONAL INPUTS:
%     Optional 'string', value pairs that may be used to modify the default
%     computational settings.
%
%     'DataOrder' - String input which sets the data order (default =
%                   'ty'). Valid inputs are 'ty' and 'yt'.
%     'Interp'    - string input controlling the interpolation method
%                   used by interp2 in the reconstruction (default =
%                   '*nearest').
%     'Plot'      - Boolean controlling whether a plot of the reconstructed
%                   estimate of the initial acoustic pressure distribution
%                   is produced (default = false). 
%     'PosCond'   - Boolean controlling whether a positivity condition is
%                   enforced on the reconstructed estimate of the initial
%                   acoustic pressure distribution (default = false).
%
% OUTPUTS:
%     p_xy        - image (indexed as x, y) 
%
% ABOUT:
%     author      - Bradley Treeby and Ben Cox
%     date        - 11th January 2009
%     last update - 21st June 2017



% avoid interp3 error in Matlab R2012b
matlab_ver = ver('matlab'); 
if strcmp(matlab_ver.Release, '(R2012b)')
    error('This function will not run in Matlab R2012b because of a bug in Matlab''s function interp2. It does work in prior and subsequent versions of Matlab.')
end

% start timer
tic;

% define defaults
num_req_inputs = 4;
data_order = 'ty';
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
                if ~strcmp(data_order, 'ty') && ~strcmp(data_order, 'yt')
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
                error('Unknown optional input.');
        end
    end
end

% reorder the data if needed (p_ty)
if strcmp(data_order, 'yt')
    p = p.';
end

% mirror the time domain data about t = 0 to allow the cosine transform to
% be computed using an FFT (p_ty)
p = [flip(p, 1); p(2:end, :)];

% extract the size of mirrored input data
[Nt, Ny] = size(p);

% update command line status
disp('Running k-Wave line reconstruction...');
disp(['  grid size: ' num2str(Ny) ' by ' num2str((Nt + 1) / 2) ' grid points']);
disp(['  interpolation mode: ' interp_method]);

% create a computational grid that is evenly spaced in w and ky, where 
% Nx = Nt and dx = dt*c
kgrid = kWaveGrid(Nt, dt * c, Ny, dy);

% from the grid for kx, create a computational grid for w using the
% relation dx = dt*c; this represents the initial sampling of p(w, ky)
w = c .* kgrid.kx;

% remap the computational grid for kx onto w using the dispersion
% relation 
% w/c = (kx^2 + ky^2)^1/2       (photoacoustics)
% w/c = (kx^2 + ky^2)/(2*kx)    (planewave ultrasound)
% This gives an w grid that is
% evenly spaced in kx. This is used for the interpolation from p(w, ky)
% to p(kx, ky). Only real w is taken to force kx (and thus x) to be
% symmetrical about 0 after the interpolation. 
%w_new = c .* kgrid.k;                      % photoacoustics
kgrid = kWaveGrid(Nt, dt*c/2, Ny, dy);      % bug fix for pwUS 15 October 2020
w_new = c .* kgrid.k.^2 ./ (2 * kgrid.kx) ; % planewave US
w_new(kgrid.kx==0) = 0;                     % planewave US

% calculate the scaling factor using the value of kx, where
% kx = sqrt( (w/c).^2 - kgrid.ky.^2 ) and then (in PA case) manually
% replacing the DC value with its limit otherwise NaN results 
% sf = c.^2 .* sqrt( (w ./ c).^2 - kgrid.ky.^2) ./ (2 .* w); % photoacoustics
% sf(w == 0 & kgrid.ky == 0) = c ./ 2;                       % photoacoustics
sf = kgrid.k + sqrt(kgrid.k.^2 - kgrid.ky.^2);               % planewave US

% compute the FFT of the input data p(t, y) to yield p(w, ky) and scale
p = sf .* fftshift(fftn(ifftshift(p)));

% remove unused variables
clear sf;

% exclude the inhomogeneous part of the wave
p(abs(w) < abs(c * kgrid.ky)) = 0;

% compute the interpolation from p(w, ky) to p(kx, ky)and then force to be
% symmetrical 
p = interp2(kgrid.ky, w, p, kgrid.ky, w_new, interp_method);

% remove unused variables
clear kgrid w;

% set values outside the interpolation range to zero
p(isnan(p)) = 0;

% compute the inverse FFT of p(kx, ky) to yield p(x, y)
p = real(fftshift(ifftn(ifftshift(p))));

% remove the left part of the mirrored data which corresponds to the
% negative part of the mirrored time data
p = p( ((Nt + 1) / 2):Nt, :);


% correct the scaling - the forward FFT is computed with a spacing of dt
% and the reverse requires a spacing of dy = dt*c, the reconstruction
% assumes that p0 is symmetrical about y, and only half the plane collects
% data (first approximation to correcting the limited view problem)
p = 2 * 2 * p ./ c;

% enfore positivity condition
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
    
    % select suitable axis scaling factor
    [~, scale, prefix] = scaleSI(max([x_axis(end), y_axis(end)])); 
    
    % select suitable plot scaling factor
    plot_scale = max(p(:));
    
    % create the figure
    figure;
    imagesc(y_axis * scale, x_axis * scale, p, [-plot_scale, plot_scale]);
    axis image;
    colormap(getColorMap);
    xlabel(['Sensor Position [' prefix 'm]']);
    ylabel(['Depth [' prefix 'm]']);
    colorbar;
    
end