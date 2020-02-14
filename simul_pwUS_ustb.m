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
pul.plot([],'broadband pulse')

% sequency generation
