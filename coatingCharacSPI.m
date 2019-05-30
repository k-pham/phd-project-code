% spin off from coatingCharac.m to assess NLprop and freqSpec from SPI time
% series taken with varying pulse energies/# of averages/different coatings

%% freqSpec and NLprop: from highly averaged SPI at different pulse energies

file_dir = '../data/coatingCharac/180221/';

file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_3.txt';

% file_name = 'ULTRA3[176us]_diffuser1b05_CNT[perspex]_allharddielectric_1.txt';
% file_name = 'ULTRA3[176us]_diffuser1b05_CNT[perspex]_allharddielectric_2.txt';
% file_name = 'ULTRA3[176us]_diffuser1b05_CNT[perspex]_allharddielectric_3.txt';

% file_name = 'ULTRA3[190us]_diffuser1b05_CNT[perspex]_allharddielectric_1.txt';
% file_name = 'ULTRA3[190us]_diffuser1b05_CNT[perspex]_allharddielectric_2.txt';
% file_name = 'ULTRA3[190us]_diffuser1b05_CNT[perspex]_allharddielectric_3.txt';

% file_name = 'ULTRA3[220us]_diffuser1b05_CNT[perspex]_allharddielectric_1.txt';
% file_name = 'ULTRA3[220us]_diffuser1b05_CNT[perspex]_allharddielectric_2.txt';
% file_name = 'ULTRA3[220us]_diffuser1b05_CNT[perspex]_allharddielectric_3.txt';

% file_name = 'ULTRA3[237us]_diffuser1b05_CNT[perspex]_allharddielectric_1.txt';
% file_name = 'ULTRA3[237us]_diffuser1b05_CNT[perspex]_allharddielectric_2.txt';
% file_name = 'ULTRA3[237us]_diffuser1b05_CNT[perspex]_allharddielectric_3.txt';

% file_name = 'ULTRA3[258us]_diffuser1b05_CNT[perspex]_allharddielectric_1.txt';
% file_name = 'ULTRA3[258us]_diffuser1b05_CNT[perspex]_allharddielectric_2.txt';
% file_name = 'ULTRA3[258us]_diffuser1b05_CNT[perspex]_allharddielectric_3.txt';

% file_name = 'ULTRA3[275us]_diffuser1b05_CNT[perspex]_allharddielectric_1.txt';
% file_name = 'ULTRA3[275us]_diffuser1b05_CNT[perspex]_allharddielectric_2.txt';
% file_name = 'ULTRA3[275us]_diffuser1b05_CNT[perspex]_allharddielectric_3.txt';

% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_avg100_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_avg100_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_avg100_3.txt';

% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-1_y-1_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-1_y-1_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-1_y-1_3.txt';

% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-2_y-1_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-2_y-1_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-2_y-1_3.txt';

% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-1_y1_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-1_y1_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_CNT[perspex]_allharddielectric_x-1_y1_3.txt';

% dt = 4e-9;
% t_min = 100;
% t_max = 300;
% freqSpecSGLsingle(file_dir,file_name,1/dt,t_min,t_max)

t_0 = 5e-6;
viewSGLsingle(file_dir,file_name,t_0)

%% freqSpec: improvement with more averages?

file_dir = '../data/coatingCharac/180227/';

file_name = 'freqSpec_CNT[perspex]_143us_avg';
array_avg = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024];
dt = 4e-9;
t_min = 100;
t_max = 300;

file_name_complete = [ file_name '1.txt' ];
freqSpecSGLsingle(file_dir,file_name_complete,1/dt,t_min,t_max)

for num_avg = array_avg(2:end)
    file_name_complete = [ file_name num2str(num_avg) '.txt' ];
    disp(file_name_complete)
    freqSpecSGLsingle(file_dir,file_name_complete,1/dt,t_min,t_max)
    pause
    drawnow
end
legend('avg1','avg2','avg4','avg8','avg16','avg32','avg64','avg128','avg256','avg512','avg1024')


%% NLprop: going to lower energies & more averages during SPI

clear all
close all

file_dir = '../data/coatingCharac/180227/';

file_name = 'NLprop_CNT[perspex]_';
% file_name = 'NLprop_SP[perspex]_';
% file_name = 'NLprop_CNT[BA59]_';
array_energy = [143, 170, 195, 215, 230, 240, 250, 260, 270, 280, 290]; %, 296];
dt = 4e-9;
t_min = 100;
t_max = 180;
t_0 = 5e-6;

for energy = array_energy(1)
    file_name_complete = [ file_name num2str(energy) 'us.txt' ];
    disp(file_name_complete)
    
    viewSGLsingle(file_dir,file_name_complete,t_0,'Norm',false)
    freqSpecSGLsingle(file_dir,file_name_complete,1/dt,t_min,t_max,'Norm',true)
    
    pause
    drawnow
end

legend('143 \mu s','170 \mu s','195 \mu s','215 \mu s','230 \mu s','240 \mu s', ...
    '250 \mu s','260 \mu s','270 \mu s','280 \mu s','290 \mu s') %,'296 \mu s')

legend('CNT on perspex','SP on perspex','CNT on sensor (BA59)')


%% pulse shape of two different Q switch delays ND filtered to match energy

file_dir = '../data/coatingCharac/180305/';

file_name = 'pulseShape_143us_ND05_CNT[perpex]_1.txt';
file_name = 'pulseShape_143us_ND05_CNT[perpex]_2.txt';
file_name = 'pulseShape_143us_ND05_CNT[perpex]_3.txt';
file_name = 'pulseShape_285us_CNT[perpex]_1.txt';
file_name = 'pulseShape_285us_CNT[perpex]_2.txt';
file_name = 'pulseShape_285us_CNT[perpex]_3.txt';

t_0 = 5e-6;
viewSGLsingle(file_dir,file_name,t_0,'Norm',false)

title('pulse shapes at different Q-switch delays ND filtered to match energy')
legend('143 \mu s + ND0.5','143 \mu s + ND0.5','143 \mu s + ND0.5', ...
    '285 \mu s','285 \mu s','285 \mu s')

% 143 us + ND0.5 has higher peak than 285 us and is thus NLprop shifted


%% freqSpec intrinsic to CNT coating (using low pressure peak with L prop)

file_dir = '../data/coatingCharac/180313/';

file_name = 'freqSpec_lowPress_270us.txt';
file_name = 'freqSpec_lowPress_280us.txt';
file_name = 'freqSpec_lowPress_290us.txt';
file_name = 'freqSpec_lowPress_297us.txt';
file_name = 'freqSpec_lowPress_300us.txt';

dt = 4e-9;
t_min = 100;
t_max = 300;
t_0 = 4e-6;

viewSGLsingle(file_dir,file_name,t_0,'Norm',true)
freqSpecSGLsingle(file_dir,file_name,1/dt,t_min,t_max,'Norm',true)

% legend('270 \mu s','280 \mu s','290 \mu s','297 \mu s','300 \mu s')


%% freqSpec intrinsic to CNT coating (using high energy pulse Nd filtered to ensure L prop)

file_dir = '../data/coatingCharac/180316/';

file_name = 'ULTRA3[143us]_diffuser1b05_ND1_CNT[perspex]_BF4_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND1_CNT[perspex]_BF4_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND1_CNT[perspex]_BF4_3.txt';
file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_BF4_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_BF4_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_BF4_3.txt';
file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_BF4_4.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_BF4_5.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_BF4_6.txt';

dt = 4e-9;
t_min = 150;
t_max = 350;
t_0 = 3.5e-6;

viewSGLsingle(file_dir,file_name,t_0,'Norm',true)
freqSpecSGLsingle(file_dir,file_name,1/dt,t_min,t_max,'Norm',true)

legend('143 \mu s + ND1','143 \mu s + ND0.5 further','143 \mu s + ND0.5 closer')


%% 180327 freqSpec intrinsic

file_dir = '../data/coatingCharac/180327/';

file_name = 'ULTRA3[143us]_diffuser1b05_ND1_CNT[perspex]_AHD_avg1024_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND1_CNT[perspex]_AHD_avg1024_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND1_CNT[perspex]_AHD_avg1024_3.txt';

% file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_AHD_avg1024_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_AHD_avg1024_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND05_CNT[perspex]_AHD_avg1024_3.txt';

% file_name = 'ULTRA3[143us]_diffuser1b05_ND03_CNT[perspex]_AHD_avg1024_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND03_CNT[perspex]_AHD_avg1024_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND03_CNT[perspex]_AHD_avg1024_3.txt';

% file_name = 'ULTRA3[143us]_diffuser1b05_ND0_CNT[perspex]_AHD_avg1024_1.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND0_CNT[perspex]_AHD_avg1024_2.txt';
% file_name = 'ULTRA3[143us]_diffuser1b05_ND0_CNT[perspex]_AHD_avg1024_3.txt';

% file_name = 'ULTRA3[143us]_diffuser1b05_closer_ND0_CNT[perspex]_AHD_avg1024.txt';

dt = 4e-9;
t_min = 150;
t_max = 350;
t_0 = 5e-6;

% viewSGLsingle(file_dir,file_name,t_0,'Norm',false)
freqSpecSGLsingle(file_dir,file_name,1/dt,t_min,t_max,'Norm',true)

% legend('143 \mu s + ND1','143 \mu s + ND0.5','143 \mu s + ND0.3','143 \mu s')


%% 190411 freqSpec with 350 MHz bandwidth scope

file_dir = 'D:\PROJECT\data\coatingCharac\190411\';

% with full power
file_name = 'ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_31s46m19h_11-04-19_avg1_savg144_raw.txt';
% file_name = 'ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_singlepoint_alignpeaks@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_17s49m19h_11-04-19_avg1_savg144_raw_jitter-corrected.txt.';

% with ND1
% accidentally using dt = 0.4 ns
% file_name = 'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[0ns]_48s09m20h_11-04-19_avg1_savg1024_raw.txt';
% file_name = 'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[0ns]_48s09m20h_11-04-19_avg1_savg1024_raw_jitter-corrected.txt';
% using dt = 0.8 ns
file_name = 'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_22s12m20h_11-04-19_avg1_savg1024_raw.txt';
% file_name = 'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_22s12m20h_11-04-19_avg1_savg1024_raw_jitter-corrected.txt';

t_0 = 0e-6;

% viewSGLsingle(file_dir,file_name,t_0,'Norm',false)
% pause

dt = 0.8e-9;
t_mins = 3000;%[3000, 3500, 3750, 3900];
t_maxs = 5000;%[5000, 4500, 4250, 4100];

% for ND1 TOA slightly later
t_mins = t_mins + 40; t_maxs = t_maxs + 40;

for idx = 1:length(t_mins)
    freqSpecSGLsingle(file_dir,file_name,1/dt,t_mins(idx),t_maxs(idx),'Norm',true)
end


%% 190507 freqSpec with 350 MHz bandwidth scope

file_dir = 'D:\PROJECT\data\coatingCharac\190507\';

% full power
% file_name = 'ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_15s41m20h_07-05-19_avg1_savg144_raw.txt';
% file_name = 'ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_15s41m20h_07-05-19_avg1_savg144_raw_jitter-corrected.txt';

% ND1 filter
% file_name = 'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_51s04m21h_07-05-19_avg1_savg1024_raw.txt';
% file_name = 'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_51s04m21h_07-05-19_avg1_savg1024_raw_jitter-corrected.txt';

% ND1 filter & distance
% file_name = 'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint_dist2@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_04s12m21h_07-05-19_avg1_savg1024_raw.txt';
% file_name = 'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint_dist2@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_04s12m21h_07-05-19_avg1_savg1024_raw_jitter-corrected.txt';

file_names = {'ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_15s41m20h_07-05-19_avg1_savg144_raw.txt'
              'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_51s04m21h_07-05-19_avg1_savg1024_raw.txt'
              'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint_dist2@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_04s12m21h_07-05-19_avg1_savg1024_raw.txt'};

t_0 = 0e-6;
dt = 0.8e-9;

t_mins = [3750 3790 6290];
t_maxs = [4250 4290 6790];
% for ND1 TOA slightly later +40
% for ND1 + dist TOA much later +2500

% viewSGLsingle(file_dir,file_name,t_0,'Norm',true,'timeAxis',false)
% pause

for idx = 2
    freqSpecSGLsingle(file_dir,file_names{idx},1/dt,t_mins(idx),t_maxs(idx),'Norm','peak','correct4PD',false)
end


%% 190517 repeat freqSpec with 350 MHz bandwidth scope with same number avg

file_dir = 'D:\PROJECT\data\coatingCharac\190517\';

file_names = {'ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_00s33m17h_17-05-19_avg1_savg1024_raw.txt'
              'ULTRA3[143us]_diffuser05b_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_00s33m17h_17-05-19_avg1_2D_raw.SGL'
              'ULTRA3[143us]_diffuser05b_ND1_CNT[perspex]_AHD1_singlepoint@0nm_t0[0]_dx[1µm]_dy[1µm]_dt[1ns]_44s38m17h_17-05-19_avg1_savg1024_raw.txt'};

t_0 = 0e-6;
dt = 0.8e-9;

t_mins = [3750 3750 3790];
t_maxs = [4250 4250 4290];
linecolours = {'m--' 'r' 'b'};

for idx = 1:3
    freqSpecSGLsingle(file_dir,file_names{idx},1/dt,t_mins(idx),t_maxs(idx), ...
        'NormTime', true, ...
        'NormFreq', true, ...
        'correct4PD', true, ...
        'LineColour',linecolours{idx} )
end

% legend('NL - uncorrected','NL - corrected for PD')
% legend('L - uncorrected','L - corrected for PD')
% legend('NL - corrected for PD','L - corrected for PD')

% legend('NL - 1024 avg','NL - 12 avg','L - 1024 avg')
