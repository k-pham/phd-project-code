% simulate simple pwUS experiment with a few scatterers using ultrasound
% toolbox - ustb

% phantom
pha = uff.phantom();
pha.sound_speed = 1500;
pha.points = [0, 0, 5e-3, 1];
fig_handle = pha.plot();

% probe (source)
prb = uff.linear_array();
prb.N = 150;
prb.pitch = 100e-6;
prb.element_width = prb.pitch;  % set width = pitch to have planar array
prb.element_height = 15e-3;
prb.plot(fig_handle);

% pulse
pul = uff.pulse();
pul.center_frequency = 1e6;
pul.fractional_bandwidth = 80;
pul.plot([],'broadband pulse');

% sequency generation
num_angle = 1;
angles = 0;
seq = uff.wave();
for n = 1 : num_angle
    seq(n) = uff.wave();
    
    seq(n).source.azimuth = angles(n);
    seq(n).source.distance = Inf;
    
    seq(n).probe = prb;
    
    seq(n).sound_speed = pha.sound_speed;
    
    % show source
    fig_handle  = seq(n).source.plot(fig_handle);
end

% simulator
sim = fresnel();          % how to use field 2 ? do i want to use field 2 ?

sim.phantom = pha;
sim.pulse = pul;
sim.probe = prb;
sim.sequence = seq;
sim.sampling_frequency = 40e6;

channel_data = sim.go();

% scan area
sca = uff.linear_scan();
sca.x_axis = linspace(-3e-3, 3e-3, 200).';
sca.z_axis = linspace(3e-3, 7e-3, 200).';
sca.plot(fig_handle, 'Scenario');

% beamformer
bmf = beamformer();
bmf.channel_data = channel_data;
bmf.scan = sca;

bmf.receive_apodization.window = uff.window.tukey50;    % is receive apodized in our case ?
bmf.receive_apodization.f_number = 1.7;                 % what is the f number ?
bmf.receive_apodization.apex.distance = Inf;

bmf.transmit_apodization.window = uff.window.tukey50;
bmf.transmit_apodization.f_number = 1.7;                % what is the f number ?
bmf.transmit_apodization.apex.distance = Inf;

% beamforming
b_data = bmf.go({ process.das_matlab() process.coherent_compounding() });

b_data.plot()




