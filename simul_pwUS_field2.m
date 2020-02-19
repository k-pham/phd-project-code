%% Computation of a pw dataset with Field II and image recon in k-space
%
% This example shows how to create a Field II simulation of Coherent Plane
% Wave Compounded (CPWC) imaging into a USTB channel_data object and
% reconstruct image with kspaceLineRecon

%% Clear old workspace and close old plots

clear all;
close all;

%% Basic Constants

c0=1500;     % Speed of sound [m/s]
fs=250e6;    % Sampling frequency [Hz]
dt=1/fs;     % Sampling step [s] 

%% field II initialisation

field_init(0);
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

%% Transducer definition L11-4v, 128-element linear array transducer

probe = uff.linear_array();
f0                      = 10e6;            % Transducer center frequency [Hz]
% lambda                  = c0/f0;           % Wavelength [m]
probe.element_height    = 15e-3;            % Height of element [m]
probe.pitch             = 0.100e-3;        % probe.pitch [m]
% kerf                    = 0.03e-03;        % gap between elements [m]
probe.element_width     = probe.pitch;     % Width of element [m]
% lens_el                 = 20e-3;           % position of the elevation focus
probe.N                 = 150;             % Number of elements
pulse_duration          = 2.5;             % pulse duration [cycles]

%% Pulse definition

pulse = uff.pulse();
pulse.fractional_bandwidth = 8;        % probe bandwidth [1]
pulse.center_frequency = f0;
t0 = (-1/pulse.fractional_bandwidth/f0): dt : (1/pulse.fractional_bandwidth/f0);


impulse_response = gauspuls(t0, f0, pulse.fractional_bandwidth);
impulse_response = impulse_response-mean(impulse_response); % To get rid of DC

te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2+1;   

% We display the pulse to check that the lag estimation is on place 
% (and that the pulse is symmetric)

figure;
plot((0:(length(two_way_ir)-1))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
plot((0:(length(two_way_ir)-1))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
legend('2-ways pulse','Envelope','Estimated lag');
title('2-ways impulse response Field II');








