clear all

file_dir = '../data/heatingEffects';

%file_name = '170707/USTplanar3M5_laserON.txt';
%file_name = '170710/USTplanar3M5_laserON_x2_y0.txt';
%file_name = '170710/USTplanar3M5_laserON_x1_y0.txt';
%file_name = '170710/USTplanar3M5_laserON_x0_y0.txt';
%file_name = '170710/USTplanar3M5_laserON_x1_y0_2.txt';
%file_name = '170710/USTplanar3M5_laserON_x1_y0_3.txt';
%file_name = '170710/USTplanar3M5_laserON_x1_y0_4.txt';

% time scale run 60.1 seconds
%file_name = '170711/timescale_500dt.txt';

% reach actual SS: leave laser ON for long then OFF for long
%file_name = '170711/USTplanar3M5_laserON_x0_y0_longSS.txt';

% measure UST in thermal equilibrium
%file_name = '170711/USTplanar3M5_laserON_x0_y0_autotrackbias.txt';

% measure laser gen US reflection in thermal equilibrium
%file_name = '170717/laserGenUS_waterReflection.txt';

% long term AC and DC triggered on laser
%file_name = '170721/laserGenUS_waterReflection_short1.txt';
%file_name = '170721/laserGenUS_waterReflection_short2.txt';
file_name = '170721/laserGenUS_waterReflection_long.txt';

%file_name = '171003/heatingEffects_BA56_longTerm_trigger[pulser].txt';
%file_name = '171003/heatingEffects_BA56_pulse2pulse_trigger[laser].txt';
%file_name = '171003/heatingEffects_BA56_pulse2pulse_trigger[laser]_2.txt';

%file_name = '171006/heatingEffects_BA56_pulse2pulse_trigger[laser].txt';
%file_name = '171006/heatingEffects_BA56_pulse2pulse_trigger[laser]_2.txt';

%file_name = '171027/heatingEffects_BA52[thin]_pulse2pulse_trigger[laser].txt';
%file_name = '171027/heatingEffects_BA52[thin]_pulse2pulse_trigger[laser]_2.txt';

%file_name = '171027/heatingEffects_BA56[thick]_pulse2pulse_trigger[laser].txt';
%file_name = '171027/heatingEffects_BA56[thick]_pulse2pulse_trigger[laser]_2.txt';

% in BEAM CENTRE !!! and changing PRF of laser:
%file_name = '171106/heatingEffects_BA56_pulse2pulse_beamcentre_20Hz.txt';
%file_name = '171106/heatingEffects_BA56_pulse2pulse_beamcentre_20Hz_2.txt';
%file_name = '171106/heatingEffects_BA56_pulse2pulse_beamcentre_10Hz.txt';
%file_name = '171106/heatingEffects_BA56_pulse2pulse_beamcentre_05Hz.txt';


file_path = [file_dir '/' file_name];

%data = importdata(filepath);

fileID = fopen(file_path,'r');
formatSpec = '%f %f %f';
sizeData = [3 Inf];
data = fscanf(fileID,formatSpec,sizeData);
fclose(fileID);
data = data';

clear file_dir and file_name and file_path

time = data(:,1);
vAC  = data(:,2);
vDC  = data(:,3);

vDC = vDC + 1;


%% hE-p2p: plot raw DC and AC consecutive acquisitions

dt_samples = 100e-9;
num_samples = 700000;
realtime = linspace(0,dt_samples*num_samples,num_samples); %in seconds
realtime = realtime - 7e-3;
realtime = realtime*1e6;    % in us

% Fig 5b heating pulse to pulse
figure
plot(realtime,vDC(1:num_samples),'r')
    xlabel('Time [\mus]')
    ylabel('Sensor signal [V]')
    xlim([-50,150])
    ylim([0.75,1.15])
    set(gca,'FontSize',14)
    set(gca,'FontName','Times New Roman')

% zoom in
figure
plot(realtime,vDC(1:num_samples),'r')
    xlabel('Time [\mus]')
    ylabel('Sensor signal [V]')
    xlim([-1,9])
    ylim([0.75,1.15])
    set(gca,'FontSize',14)
    set(gca,'FontName','Times New Roman')

