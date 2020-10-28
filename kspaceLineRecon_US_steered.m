function F = kspaceLineRecon_US_steered(p, dy, dt, c, a, varargin)
%KSPACELINERECON_US_STEERED 2D linear FFT reconstruction for plane wave
%                   ultrasound reflection imaging with steering angle
%
% DESCRIPTION
%     based on kspaceLineRecon_US
% 
% USAGE:
%     F_xy = kspaceLineRecon_US_steered(p_ty, dy, dt, c, angle)
%     F_xy = kspaceLineRecon_US_steered(p_ty, dy, dt, c, angle, ...)
%
% INPUTS:
%     p_ty        - pressure time-series recorded over an evenly spaced
%                   array of sensor points on a line (indexed as t, y)
%     dy          - spatial step [m]
%     dt          - time step [s]
%     c           - acoustically-homogeneous sound speed [m/s]
%     a           - steering angle of transmitted plane wave [rad]
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
%     F_xy        - image of object function (indexed as x, y) 
%
% ABOUT:
%     author      - Bradley Treeby and Ben Cox
%     date        - 11th January 2009
%     update      - 21st June 2017 (Ben - mod for pwUS)
%     update      - 16 October 2020 (Khoa - corrections for pwUS)
%     update      - 23 October 2020 (Khoa - rename output & tidy comments &
%                                           explicitly distinquish between
%                                           receive and object k space)
%     last update - 27 October 2020 (Khoa - mod for steered pwUS)



% avoid interp3 error in Matlab R2012b
matlab_ver = ver('matlab'); 
if strcmp(matlab_ver.Release, '(R2012b)')
    error('This function will not run in Matlab R2012b because of a bug in Matlab''s function interp2. It does work in prior and subsequent versions of Matlab.')
end

% start timer
tic;

% define defaults
num_req_inputs = 5;
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

% create a computational grid that is evenly spaced in w and ky (the
% receive k-space), where the kx component is w/c, so Nx = Nt and dx = dt*c
kgrid_rec = kWaveGrid(Nt, dt*c, Ny, dy);

% from the grid for kx, create a computational grid for w using the
% relation dx = dt*c; this represents the initial sampling of p(w, ky)
w = c .* kgrid_rec.kx;

% calculate the scaling factor using the value of kz, where
% kz = sqrt( (w/c)^2 - ky^2 )
% sf = c.^2 .* sqrt( (w ./ c).^2 - kgrid.ky.^2) ./ (2 .* w); % photoacoustics
% sf(w == 0 & kgrid.ky == 0) = c ./ 2;                       % photoacoustics
% sf = kgrid.k + sqrt(kgrid.k.^2 - kgrid.ky.^2);             % planewave US
sf = sqrt( (w./c).^2 - kgrid_rec.ky.^2 );                    % correction for pwUS 16 October 2020

% compute the FFT of the input data p(t, y) to yield p(w, ky) and scale
% with sf to give F(w, ky)
F = sf .* fftshift(fftn(ifftshift(p)));

% remove unused variables
clear sf;

% exclude the inhomogeneous part of the wave
F(abs(w) < abs(c * kgrid_rec.ky)) = 0;

% create a new computational grid that is evenly spaced in kz' and ky' (the
% object k-space) where factor of 1/2 in dz' due to reflection imaging
% (this step is not necessary for photoacoustics, since object kgrid
% happens to be the same as receive kgrid)
kgrid_obj = kWaveGrid(Nt, dt*c/2, Ny, dy);      % bug fix for pwUS 15 October 2020

% remap the computational grid for kz' onto w using the dispersion relation
% w/c = (ky'^2 + kz'^2)^1/2                            (photoacoustics)
% w/c = (ky'^2 + kz'^2)/(2*kz')                        (planewave ultrasound)
% w/c = (ky'^2 + kz'^2)/(2*ky'*sin(a)+2*kz'*cos(a))    (steered pw US)
% ky  = ky' - w*sin(a)/c;                              (steered pw US)
% This gives an w grid that is evenly spaced in kz'. This is used for the 
% interpolation from F(w, ky) to F(kz', ky'). Only real w is taken to force
% kz' (and thus z) to be symmetrical about 0 after the interpolation.
% w_new = c .* kgrid.k;                               % photoacoustics
% w_new = c .* kgrid_obj.k.^2 ./ (2 * kgrid_obj.kx) ; % planewave US
% w_new(kgrid_obj.kx==0) = 0;                         % planewave US
denominator = 2*kgrid_obj.ky.*sin(a) + 2*kgrid_obj.kx.*cos(a);   % steered pw US
w_new = c .* kgrid_obj.k.^2 ./ denominator;                      % steered pw US
w_new(denominator==0) = 0;                                       % steered pw US
ky_new = kgrid_obj.ky - w_new .* sin(a) ./ c;                    % steered pw US

% compute the interpolation from F(w, ky) to F(kz', ky') and then force to
% be symmetrical
F = interp2(kgrid_rec.ky, w, F, ky_new, w_new, interp_method);

% remove unused variables
clear kgrid_rec kgrid_obj w;

% set values outside the interpolation range to zero
F(isnan(F)) = 0;

% compute the inverse FFT of p(kx, ky) to yield p(x, y)
F = real(fftshift(ifftn(ifftshift(F))));

% remove the left part of the mirrored data which corresponds to the
% negative part of the mirrored time data
F = F( ((Nt + 1) / 2):Nt, :);

% correct the scaling - the forward FFT is computed with a spacing of dt
% and the reverse requires a spacing of dy = dt*c, the reconstruction
% assumes that p0 is symmetrical about y, and only half the plane collects
% data (first approximation to correcting the limited view problem)
F = 2 * 2 * F ./ c;

% enforce positivity condition
if positivity_cond
    disp('  applying positivity condition...');
    F(F < 0) = 0;
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
    plot_scale = max(F(:));
    
    % create the figure
    figure;
    imagesc(y_axis * scale, x_axis * scale, F, [-plot_scale, plot_scale]);
    axis image;
    colormap(getColorMap);
    xlabel(['Sensor Position [' prefix 'm]']);
    ylabel(['Depth [' prefix 'm]']);
    colorbar;
end

end