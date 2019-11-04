% USimagingPhantoms.m:      parameters and file locations
% USimagingRecon.m:         run reconstruction, save .mat files, plot
% USimagingAnalysis.m:      image analysis


% dim                   dimension of imaging: 2=line scan, 3=plane scan
% dx                    grid point spacing in the x direction [m]
% dy                    grid point spacing in the y direction [m]
% trigger_delay         delay of acquisition post-trigger [s]
% samples_cut_off       number of samples to zero at the beginning
% samples_t0_correct    number of samples to add to correct for t0
% c0                    sound speed [m/s]
% file_name             'date\phantomID@[...].SGL'

file_dir = '..\data\imagingUS\';
% file_dir = '..\data\imagingDM\';


%% 170724 resolution80um (glass fibres in water)

% dim = 2;
% dx = 0.1e-3;
% trigger_delay = 0;

% samples_cut_off = 50;
% samples_t0_correct = 0;
% c0 = 1484;

% file_name = '170724/resolution80um_waterReflection@0nm_t0[-1500]_dx[20µm]_dy[0µm]_dt[4ns]_12s16m19h_24-07-17_avg1_1D_raw.SGL';
% file_name = '170724/resolution80um@0nm_t0[-1750]_dx[20µm]_dy[0µm]_dt[4ns]_43s22m19h_24-07-17_avg1_1D_raw.SGL';
% file_name = '170724/resolution80um_short@0nm_t0[-1750]_dx[20µm]_dy[0µm]_dt[4ns]_53s24m19h_24-07-17_avg1_1D_raw.SGL';


%% 170801 resolution27um (tungsten wires in water)

% dim = 2;
% dx = 0.02e-3;
% dt = 4e-9;
% trigger_delay = 0;
% 
% samples_cut_off = 50;
% samples_t0_correct = -9;
% c0 = 1479;
% 
% % file_name = '170801/resolution27_waterReflection@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_19s48m13h_01-08-17_avg16_1D_raw.SGL';
% % file_name = '170801/resolution27@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_05s22m14h_01-08-17_avg1_1D_raw.SGL';
% file_name = '170801/resolution27@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_03s37m14h_01-08-17_avg16_1D_raw.SGL';


%% 170802 layerGelWax (gel wax layers coupled with coupling gel)
 
% dim = 2;
% dx = 0.1e-3;
% trigger_delay = 0;

% samples_cut_off = 50;
% samples_t0_correct = 2;
% c0 = 1580;

% file_name = '170802/layerGelWax_5pp_2pp_5SGS@0nm_t0[0]_dx[0µm]_dy[100µm]_dt[4ns]_28s48m12h_02-08-17_avg1_1D_raw.SGL';
% file_name = '170802/layerGelWax_5pp_2pp_5SGS@0nm_t0[0]_dx[0µm]_dy[100µm]_dt[4ns]_42s50m12h_02-08-17_avg8_1D_raw.SGL';
% file_name = '170802/layerGelWax_5SGS_2pp_5pp@0nm_t0[0]_dx[0µm]_dy[100µm]_dt[4ns]_39s21m14h_02-08-17_avg8_1D_raw.SGL';


%% 170807 polymerLeaf (polymer leaf skeleton in water)

% dim = 3;
% dx = 0.1e-3;
% dy = 0.1e-3;
% dt = 4e-9;
% trigger_delay = 8e-6;
% 
% samples_cut_off = 0;
% samples_t0_correct = 2;
% c0 = 1484;
% 
% file_name = '170807/polymerLeaf_long@0nm_t0[-2000]_dx[100µm]_dy[100µm]_dt[4ns]_39s08m17h_07-08-17_avg1_2D_raw.SGL';


%% 170815 polymerLeaf (polymer leaf skeleton in water)

% dim = 3;
% dx = 0.1e-3;
% dy = 0.1e-3;
% dt = 4e-9;
% trigger_delay = 8e-6;
% 
% samples_cut_off = 0;
% samples_t0_correct = -9;
% c0 = 1484;
% 
% file_name = '170815/polymerLeaf_long4@0nm_t0[-2000]_dx[100µm]_dy[100µm]_dt[4ns]_02s48m12h_15-08-17_avg1_2D_raw.SGL';


%% 180416 resolution27um (tungsten wires in water)

% dim = 2;
% dx = 0.02e-3;
% trigger_delay = 0;

% samples_cut_off = 50;
% samples_t0_correct = 2;
% c0 = 1467;

% file_name = '180416\resolution27um_BA59[CNT]_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_14s33m17h_16-04-18_avg1_1D_raw.SGL';
% file_name = '180416\resolution27um_BA59[CNT]_avg16@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_05s47m17h_16-04-18_avg16_1D_raw.SGL';


%% 180509 resolution6um (carbon fibres in water)

% dim = 2;
% dx = 0.02e-3;
% trigger_delay = 0;
% 
% samples_cut_off = 50;
% samples_t0_correct = 2;
% c0 = 1484;
% 
% file_name = '180509\resolution6um_BA21[CNT]_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_00s41m14h_09-05-18_avg1_1D_raw.SGL';
% file_name = '180509\resolution6um_BA21[CNT]_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_39s48m14h_09-05-18_avg10_1D_raw.SGL';


%% 180515 resolution6um (carbon fibres in water) move around in space

% dim = 2;
% dx = 0.02e-3;
% trigger_delay = 0;
% 
% samples_cut_off = 100;
% samples_t0_correct = -20;
% c0 = 1484;
% 
% file_names = {
%     '180515\resolution6um_BA21[CNT]_z0.5_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_13s15m21h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z0.5_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_04s24m21h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z1_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_20s32m21h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z1_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_35s40m21h_15-05-18_avg10_1D_raw.SGL';
%     '180515\resolution6um_BA21[CNT]_z1.5_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_44s52m21h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z1.5_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_39s00m22h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z2_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_33s08m22h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z2_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_04s16m22h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z2.5_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_15s26m22h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z2.5_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_48s33m22h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z3_x30.5_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_21s41m22h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z3_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_49s48m22h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z4_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_25s58m22h_15-05-18_avg10_1D_raw.SGL'
%     '180515\resolution6um_BA21[CNT]_z5_x31_avg10@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_27s08m23h_15-05-18_avg10_1D_raw.SGL'
%     };
% 
% % need to reshape to be able to loop through
% file_names = reshape(file_names,[1 length(file_names)]);
% 
% % moved to separate .m file
% % also facilitated within USimagingRecon.m


%% 180518 resolution27um (tungsten wires in water) compare BA59 vs BA21

% dim = 2;
% dx = 0.02e-3;
% trigger_delay = 0;
% 
% samples_cut_off = 100;
% samples_t0_correct = -20;
% c0 = 1484;
% 
% file_names = {
% %     '180518\resolution27um_BA21[CNT]_z1_x28_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_37s38m18h_18-05-18_avg1_1D_raw.SGL'
% %     '180518\resolution27um_BA21[CNT]_z1_x28_avg2@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_38s41m18h_18-05-18_avg2_1D_raw.SGL'
% %     '180518\resolution27um_BA21[CNT]_z1_x28_avg4@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_13s46m18h_18-05-18_avg4_1D_raw.SGL'
%     '180518\resolution27um_BA21[CNT]_z1_x28_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_44s32m19h_18-05-18_avg1_1D_raw.SGL'
%     '180518\resolution27um_BA59[CNT]_z1_x28_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_27s07m22h_18-05-18_avg1_1D_raw.SGL'
% %     '180518\resolution27um_BA59[CNT]_z1_x28_avg2@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_45s39m22h_18-05-18_avg2_1D_raw.SGL'
% %     '180518\resolution27um_BA59[CNT]_z1_x28_avg4@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_55s52m22h_18-05-18_avg4_1D_raw.SGL'
%     };
% 
% % need to reshape to be able to loop through
% file_names = reshape(file_names,[1 length(file_names)]);


%% 180519 resolution27um (tungsten wires fibres in water) 

% dim = 2;
% dx = 0.02e-3;
% trigger_delay = 0;
% 
% samples_cut_off = 100;
% samples_t0_correct = -16;
% c0 = 1484;
% 
% file_names = {
%     '180519\resolution27um_BA59[CNT]_z1_x28_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_49s19m17h_19-05-18_avg1_1D_raw.SGL'
%     '180519\resolution27um_BA59[CNT]_z1_x28_avg2@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_00s23m17h_19-05-18_avg2_1D_raw.SGL'
%     '180519\resolution27um_BA59[CNT]_z1_x28_avg4@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_45s28m17h_19-05-18_avg4_1D_raw.SGL'
%     '180519\resolution27um_BA59[CNT]_z1_x28_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_13s46m17h_19-05-18_avg1_1D_raw.SGL'
%     '180519\resolution27um_BA59[CNT]_z1_x28_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_11s10m18h_19-05-18_avg1_1D_raw.SGL'
%     '180519\resolution27um_BA59[CNT]_z1.5_x28_avg1@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_42s21m18h_19-05-18_avg1_1D_raw.SGL'
%     };
% 
% % need to reshape to be able to loop through
% file_names = reshape(file_names,[1 length(file_names)]);


%% 180522 deadMouse in water

% dim = 3;
% dx = 0.1e-3;
% dy = 0.1e-3;
% dt = 4e-9;
% trigger_delay = 0;
% 
% samples_cut_off = 100;
% samples_t0_correct = -20;
% c0 = 1540;
% 
% % file_name = '180522\deadMouse_BA21[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_33s57m17h_22-05-18_avg1_2D_raw.SGL';
% file_name = '180522\deadMouse_on_side_BA21[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_18s07m19h_22-05-18_avg1_2D_raw.SGL';


%% 180524 resolution27um (tungsten wires in water) moved in space

% dim = 2;
% dx = 0.02e-3;
% % trigger_delay = 11.5e-6;
% 
% samples_cut_off = 100;
% samples_t0_correct = -20;
% c0 = 1484; % 1486
% 
% diff_line_scans = {
% %     '180524\resolution27um_BA21[CNT]_z0.5_x25@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_33s22m15h_24-05-18_avg1_1D_raw.SGL' 0
% %     '180524\resolution27um_BA21[CNT]_z0.5_x25@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_55s30m15h_24-05-18_avg1_1D_raw.SGL' 0
% %     '180524\resolution27um_BA21[CNT]_z0.5_x25@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_08s38m15h_24-05-18_avg1_1D_raw.SGL' 0
% %     '180524\resolution27um_BA21[CNT]_z0.5_x24@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_04s44m15h_24-05-18_avg1_1D_raw.SGL' 0
% %     '180524\resolution27um_BA21[CNT]_z0.5_x26@0nm_t0[0]_dx[20µm]_dy[0µm]_dt[4ns]_19s47m15h_24-05-18_avg1_1D_raw.SGL' 0
% %     '180524\resolution27um_BA21[CNT]_z3.5_x26@0nm_t0[-1000]_dx[20µm]_dy[0µm]_dt[4ns]_48s56m15h_24-05-18_avg1_1D_raw.SGL' 4e-6
% %     '180524\resolution27um_BA21[CNT]_z3_x26@0nm_t0[-875]_dx[20µm]_dy[0µm]_dt[4ns]_06s05m16h_24-05-18_avg1_1D_raw.SGL' 3.5e-6
% %     '180524\resolution27um_BA21[CNT]_z3_x25@0nm_t0[-875]_dx[20µm]_dy[0µm]_dt[4ns]_16s10m16h_24-05-18_avg1_1D_raw.SGL' 3.5e-6
% %     '180524\resolution27um_BA21[CNT]_z3_x24@0nm_t0[-875]_dx[20µm]_dy[0µm]_dt[4ns]_11s12m16h_24-05-18_avg1_1D_raw.SGL' 3.5e-6
% %     '180524\resolution27um_BA21[CNT]_z6_x24@0nm_t0[-1875]_dx[20µm]_dy[0µm]_dt[4ns]_53s18m16h_24-05-18_avg1_1D_raw.SGL' 7.5e-6
% %     '180524\resolution27um_BA21[CNT]_z6_x25@0nm_t0[-1875]_dx[20µm]_dy[0µm]_dt[4ns]_32s24m16h_24-05-18_avg1_1D_raw.SGL' 7.5e-6
% %     '180524\resolution27um_BA21[CNT]_z6_x26@0nm_t0[-1875]_dx[20µm]_dy[0µm]_dt[4ns]_28s26m16h_24-05-18_avg1_1D_raw.SGL' 7.5e-6
% %     '180524\resolution27um_BA21[CNT]_z9_x26@0nm_t0[-2875]_dx[20µm]_dy[0µm]_dt[4ns]_28s36m16h_24-05-18_avg1_1D_raw.SGL' 11.5e-6
%     '180524\resolution27um_BA21[CNT]_z9_x25@0nm_t0[-2875]_dx[20µm]_dy[0µm]_dt[4ns]_24s40m16h_24-05-18_avg1_1D_raw.SGL' 11.5e-6
% %     '180524\resolution27um_BA21[CNT]_z9_x24@0nm_t0[-2875]_dx[20µm]_dy[0µm]_dt[4ns]_09s42m16h_24-05-18_avg1_1D_raw.SGL' 11.5e-6
%     };
% 
% file_names      = diff_line_scans(:,1);
% trigger_delays  = diff_line_scans(:,2);
% 
% % need to reshape to be able to loop through
% file_names     = reshape(file_names,[1 length(file_names)]);
% trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);


%% 180524 mouseKidney (ex vivo in water)

% dim = 3;
% dx = 0.1e-3;
% dy = 0.1e-3;
% dt = 4e-9;
% trigger_delay = 0;
% 
% samples_cut_off = 100;
% samples_t0_correct = -20;
% c0 = 1560;
% 
% file_name = '180524\mouseKidney_BA21[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_44s42m17h_24-05-18_avg1_2D_raw.SGL';


%% 180525 mouseKidney (ex vivo in water covered by thin polymer film to prevent floating)

% dim = 3;
% dx = 0.1e-3;
% dy = 0.1e-3;
% dt = 20e-9;
% trigger_delay = 0;
% 
% samples_cut_off = 100;
% samples_t0_correct = -4;
% c0 = 1560;
% 
% file_name = '180525\mouseKidney_BA21[CNT]_long@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_41s24m18h_25-05-18_avg1_2D_raw.SGL';


%% 180529 optifibreKnot in water

% dim = 3;
% dx = 0.1e-3;
% dy = 0.1e-3;
% dt = 8e-9;
% trigger_delay = 8e-6;
% 
% samples_cut_off = 0;
% samples_t0_correct = -20;
% c0 = 1484;
% 
% file_name = '180529\optifibreKnot_BA21[CNT]@0nm_t0[-1000]_dx[100µm]_dy[100µm]_dt[8ns]_35s53m21h_29-05-18_avg1_2D_raw.SGL';


%% 180529 PVAphantom1 in water

% dim = 3;
% dx = 0.1e-3;
% dy = 0.1e-3;
% dt = 8e-9;
% trigger_delay = 8e-6;
% 
% samples_cut_off = 0;
% samples_t0_correct = -20;
% c0 = 1520; % ??
% 
% file_name = '180529\PVAphantom1_BA21[CNT]@0nm_t0[-1000]_dx[100µm]_dy[100µm]_dt[8ns]_25s47m22h_29-05-18_avg1_2D_raw.SGL';


%% 180531 layers of PVA coupled by water (trolley scanner)

% dim = 3;
% dx = 0.1e-3;
% dy = 0.1e-3;
% dt = 5e-9;
% trigger_delay = 0;
% 
% samples_cut_off = 0;
% samples_t0_correct = -17;
% c0 = 1520; % ??
% 
% file_name = '180531\pvaLayers_BA21[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[5ns]_04s34m21h_31-05-18_avg1_2D_raw.SGL';


%% 180619 optifibreKnot_angled & polymerLeaf in water (trolley scanner - new PD board)

% file_name = '180619/optifibreKnot_angled_BA21[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_31s22m17h_19-06-18_avg1_2D_raw.SGL';
% file_name = '180619/polymerLeaf_BA21[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_05s41m19h_19-06-18_avg1_2D_raw.SGL';
% 
% [sensor_data, params] = loadSGL([file_dir file_name]);
% 
% dim = 3;
% dx = params.dx;
% dy = params.dy;
% dt = params.dt;
% trigger_delay = 8e-6;
% 
% samples_cut_off = 0;
% samples_t0_correct = -6;
% c0 = 1484;


%% 180621 polymerLeafFlat in water (moved trolley and ultra2 into trung's lab, scrambled fibre)

% file_name = '180621\polymerLeafFlat_BA21[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_46s45m19h_21-06-18_avg1_2D_raw.SGL';
% 
% dim = 3;
% 
% trigger_delay = 8e-6;
% 
% samples_cut_off = 0;
% samples_t0_correct = -6;
% c0 = 1484;


%% 180625 polymerLeafFlat in water - good sensor BK31 (but heatingEffects?)

% file_name = '180625\polymerLeafFlat_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_08s15m19h_25-06-18_avg1_2D_raw.SGL';
% 
% dim = 3;
% 
% trigger_delay = 8e-6;
% 
% samples_cut_off = 0;
% samples_t0_correct = -6;
% c0 = 1484;


%% 180626 3D phantoms - good sensor BK31
% optifibreKnot mildly angled
% repeat with less angle/more central/smaller knot
% polymerLeaf in water
% polymerLeafFlat in water
% polymerHairBundle in water

% file_name = '180626\optifibreKnot_angled_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_56s21m13h_26-06-18_avg1_2D_raw.SGL';
% file_name = '180626\optifibreKnot_angled2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_31s07m15h_26-06-18_avg1_2D_raw.SGL';
% file_name = '180626\optifibreKnot_angled3_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_44s39m16h_26-06-18_avg1_2D_raw.SGL';
% file_name = '180626\optifibreKnot_angled4_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_10s39m17h_26-06-18_avg1_2D_raw.SGL';
% file_name = '180626\optifibreKnot_angled5_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_30s17m18h_26-06-18_avg1_2D_raw.SGL';
% trigger_delay = 5e-6;

% file_name = '180626\polymerLeaf_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_39s38m19h_26-06-18_avg1_2D_raw.SGL';
% trigger_delay = 4.5e-6;

% file_name = '180626\polymerLeaf2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_51s28m20h_26-06-18_avg1_2D_raw.SGL';
% file_name = '180626\polymerLeafFlat2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_07s00m22h_26-06-18_avg1_2D_raw.SGL';
% trigger_delay = 4e-6;

% file_name = '180626\polymerHairBundle_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_16s50m22h_26-06-18_avg1_2D_raw.SGL';
% trigger_delay = 0;

% dim = 3;

% samples_cut_off = 0;
% samples_t0_correct = -6;
% c0 = 1484;


%% 180626 pva block + needles (US looking phantom)

% file_name = '180626\pvaBlockNeedles_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_55s37m23h_26-06-18_avg1_2D_raw.SGL';
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 10;
% samples_t0_correct = -6;
% c0 = 1540;


%% 180627 gel wax block + cylinders of diff scat (US looking phantom)

% file_name = '180627\gelwaxBlockCylinders1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_24s39m13h_27-06-18_avg1_2D_raw.SGL';
% trigger_delay = 0;
% samples_cut_off = 10;

% file_name = '180627\gelwaxBlockCylinders2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_05s11m14h_27-06-18_avg1_2D_raw.SGL';
% trigger_delay = 5e-6;
% samples_cut_off = 0;
% 
% dim = 3;
% 
% samples_t0_correct = -6;
% c0 = 1460;


%% 180627 turkey breast with hairs threaded through

% file_name = '180627\turkeyBreastHairs_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_31s53m14h_27-06-18_avg1_2D_raw.SGL';
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 10;
% samples_t0_correct = -6;
% c0 = 1560;


%% 180629 gel wax layers coupled by coupling gel

% file_name = '180629\gelwaxLayers_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[10ns]_39s51m16h_29-06-18_avg1_2D_raw.SGL';
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 10;
% samples_t0_correct = -6;
% c0 = 1470;
% 
% % clinical scanner:
% file_name_cli = '180629\clinical scanner\29-06-2018_17-41-25 DICOM +2 higher gain\17-38-33.dcm';

%% 180629 pork belly & turkey breast w hairs in diagonal

% file_name = '180629\porkBelly_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_12s50m21h_29-06-18_avg1_2D_raw.SGL';
% c0 = 1460; % ??

% file_name = '180629\turkeyBreastHairs2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_03s33m22h_29-06-18_avg1_2D_raw.SGL';
% c0 = 1580;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 10;
% samples_t0_correct = -6;


%% 180814 pork belly 2

% file_name = '180814\porkBelly2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_23s36m17h_14-08-18_avg1_2D_raw.SGL';
% c0 = 1460; % ??
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 10;
% samples_t0_correct = -6;


%% 180821 chicken breast with water-filled glass capillaries

% file_name = '180821\chickenBreastCapillaries_BK31@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_21s54m18h_21-08-18_avg1_2D_raw.SGL';
% c0 = 1560;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 10;
% samples_t0_correct = -5;


%% 180822 tofu block w/ holes (punched with big needle)

% file_name = '180822\tofuBlockHoles_BK31@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_09s30m15h_22-08-18_avg1_2D_raw.SGL';
% c0 = 1485;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 10;
% samples_t0_correct = -5;


%% 180828 pork belly 3

% file_name = '180828\porkBelly3_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[20ns]_26s20m19h_28-08-18_avg1_2D_raw.SGL';
% c0 = 1460;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 10;
% samples_t0_correct = -4;


%% 180919 polymerLeafFlat BF26 with AuNP coating

% file_name = '180919\polymerLeafFlat_BF26[AuNP]_600nm@0nm_t0[0]_dx[150µm]_dy[150µm]_dt[4ns]_11s59m12h_19-09-18_avg1_2D_raw.SGL';
% c0 = 1502;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -27;

% PA data for US recon (to see crosstalk)
% file_name = '180919\polymerLeafFlat_BF26[AuNP]_900nm@0nm_t0[-500]_dx[150µm]_dy[150µm]_dt[4ns]_57s13m12h_19-09-18_avg1_2D_raw.SGL';
% c0 = 1502;
% 
% trigger_delay = 2e-6;
% 
% dim = 3;
% 
% samples_cut_off = 0;
% samples_t0_correct = -27;


%% 181010 hair BF26[AuNP]

% file_name = '181010\hair_BF26[AuNP]_600nm@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_56s04m13h_10-10-18_avg1_2D_raw.SGL';
% c0 = 1502;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -28;


%% 181015 mouseKidneyInAgar BK31[CNT]

% file_name = '181015\mouseKidneyInAgar_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_19s18m18h_15-10-18_avg1_2D_raw.SGL';
% c0 = 1500;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -9;


%% 181019 mouseKidneyInAgar2/3/_flipped BK31[CNT]

% file_name = '181019\mouseKidneyInAgar2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_56s50m17h_19-10-18_avg1_2D_raw.SGL';
% file_name = '181019\mouseKidneyInAgar2_flipped_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_49s33m18h_19-10-18_avg1_2D_raw.SGL';
% file_name = '181019\mouseKidneyInAgar3_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_21s05m19h_19-10-18_avg1_2D_raw.SGL';
% file_name = '181019\mouseKidneyInAgar3_flipped_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_47s32m19h_19-10-18_avg1_2D_raw.SGL';
% c0 = 1500;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -9;


%% 181030 mouseKidneyInOptiLube BK31[CNT] - BAD IMAGE LOTS OF BUBBLES/HIGH INTENSITY NOISE

% file_name = '181030\mouseKidneyInOptiLube_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_09s27m15h_30-10-18_avg1_2D_raw.SGL';
% c0 = 1520;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -9; % try -6 -> -9


%% 181102 mouseKidneyInAgar #4 (flipped) BK31[CNT]

% file_name = '181102\mouseKidneyInAgar4_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_25s45m17h_02-11-18_avg1_2D_raw.SGL';
% file_name = '181102\mouseKidneyInAgar42_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_29s05m19h_02-11-18_avg1_2D_raw.SGL';
% file_name = '181102\mouseKidneyInAgar42_BK31[CNT]_avg5@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_30s48m20h_02-11-18_avg5_2D_raw.SGL';
% c0 = 1520;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -9; % try -6 -> -9


%% 181119 polymerLeaf (curved, flipped) BK31[CNT]

% file_name = '181119\polymerLeaf_flipped_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_28s20m12h_19-11-18_avg1_2D_raw.SGL';
% c0 = 1484;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -9;


%% 181119 hair - multiple times for averaging (SNR increased ?) BK31[CNT]

% file_name = '181119\hair_BK31[CNT]_1@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_14s28m16h_19-11-18_avg1_2D_raw.SGL';
% file_name = '181119\hair_BK31[CNT]_2@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_46s42m16h_19-11-18_avg1_2D_raw.SGL';
% file_name = '181119\hair_BK31[CNT]_3@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_36s53m16h_19-11-18_avg1_2D_raw.SGL';
% file_name = '181119\hair_BK31[CNT]_4@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_50s03m17h_19-11-18_avg1_2D_raw.SGL';
% file_name = '181119\hair_BK31[CNT]_5@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_07s14m17h_19-11-18_avg1_2D_raw.SGL';
% file_name = '181119\hair_BK31[CNT]_6@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_17s24m17h_19-11-18_avg1_2D_raw.SGL';
% c0 = 1484;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -9;


%% 181126 hair bundle - multi/16-beam by Trung BF26(?)[CNT] - 850 nm 100 Hz OPO

% file_name = '/181126 (Trung multi-beam)/Full_scan1_hair1@850nm_t0[0]_dx[106µm]_dy[106µm]_dt[17ns]_59s15m15h_26-11-18_avg1_2D_raw.SGL';
% % file_name = '/181126 (Trung multi-beam)/Full_scan1_hair2@850nm_t0[0]_dx[106µm]_dy[106µm]_dt[17ns]_38s00m16h_26-11-18_avg1_2D_raw.SGL';
% 
% c0 = 1484;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -15;


%% 181128 lymph node 1 - BK31[CNT]

% file_name = '181128/lymphNode1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_00s13m14h_28-11-18_avg1_2D_raw.SGL';
% 
% c0 = 1520;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -6;


%% 181204 agar-based tissue mimicking phantom (atmm) with orgasol - BK31[CNT]

% file_name = '181204/atmm_orgasol1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_03s08m21h_04-12-18_avg1_2D_raw.SGL';
% 
% c0 = 1544;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -6;


%% 181205 bacon layers & hair - BK31[CNT]

% file_name = '181205/baconLayersHair_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_05s14m12h_05-12-18_avg1_2D_raw.SGL';
% 
% c0 = 1540;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 0;
% samples_t0_correct = -6;


%% 181211 atmm with SiC - BK31[CNT]

% file_name = '181211/atmm_SiC1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_07s05m23h_11-12-18_avg1_2D_raw.SGL'; % moved while imaging
% file_name = '181211/atmm_SiC1_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_54s40m23h_11-12-18_avg1_2D_raw.SGL';
% 
% c0 = 1544;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -6;


%% 181214 atmm sith SiC #2 - BK31[CNT]

% file_name = '181214/atmm_SiC2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_26s27m19h_14-12-18_avg1_2D_raw.SGL';
% 
% c0 = 1544;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -6;


%% 190114 lymph node (L3) 2nd try - BK31[CNT]

% file_name = '190114/lymphNode2_BK31[CNT]@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_13s51m16h_14-01-19_avg1_2D_raw.SGL';
% 
% c0 = 1520;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -6;


%% 190122 DM polymer leaf - BF16[AuNP]

% file_name = '190122/polymerLeaf_BF16[AuNP]_600nm@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_52s22m21h_22-01-19_avg1_2D_raw.SGL';
% file_name = '190122/polymerLeaf_BF16[AuNP]_900nm_02@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_15s25m22h_22-01-19_avg1_2D_raw.SGL';
% 
% c0 = 1484;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 100;
% samples_t0_correct = -10;


%% 190124 DM atmm + india ink tubes - BF16[AuNP]

% file_name = '190124/atmm+indiaInk_BF16[AuNP]_600nm@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_35s46m19h_24-01-19_avg1_2D_raw.SGL';
% file_name = '190124/atmm+indiaInk2_BF16[AuNP]_600nm@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_35s21m21h_24-01-19_avg1_2D_raw.SGL';
% 
% c0 = 1544;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -5;


%% 190129 DM agar + ink/water tubes(PMMA) - BF16[AuNP]

% file_name = '190129/agar+inkORwater_BF16[AuNP]_600nm@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_41s54m21h_29-01-19_avg1_2D_raw.SGL';
% 
% c0 = 1544;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -5;


%% 190129 ink tubes(soft rubber) - BF16[AuNP]

% file_name = '190129/inkTubes_BF16[AuNP]_600nm@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[8ns]_31s20m23h_29-01-19_avg1_2D_raw.SGL';
% 
% c0 = 1484;
% 
% trigger_delay = 0;
% 
% dim = 3;
% 
% samples_cut_off = 50;
% samples_t0_correct = -5;


%% 190606 resolution27umPlanar - BK31[CNT]

% dim = 3;

% file_name = '190606/resolution27umPlanar_BK31[CNT]_2D@0nm_t0[0]_dx[200µm]_dy[200µm]_dt[4ns]_41s47m15h_06-06-19_avg1_2D_raw.SGL';
% file_name = '190606/resolution27umPlanar_BK31[CNT]_2D_2@0nm_t0[0]_dx[200µm]_dy[200µm]_dt[4ns]_13s05m16h_06-06-19_avg1_2D_raw.SGL';
% file_name = '190606/resolution27umPlanar_BK31[CNT]_2D_3@0nm_t0[0]_dx[200µm]_dy[200µm]_dt[4ns]_06s21m16h_06-06-19_avg1_2D_raw.SGL';
% file_name = '190606/resolution27umPlanar_BK31[CNT]_2D_4@0nm_t0[0]_dx[200µm]_dy[200µm]_dt[4ns]_12s32m16h_06-06-19_avg1_2D_raw.SGL';

% trigger_delay = 0;

% dim = 2;
% 
% diff_line_scans = {
%     '190606/resolution27umPlanar_BK31[CNT]_1@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_27s41m16h_06-06-19_avg1_1D_raw.SGL' 0
%     '190606/resolution27umPlanar_BK31[CNT]_2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_32s43m16h_06-06-19_avg1_1D_raw.SGL' 0
%     '190606/resolution27umPlanar_BK31[CNT]_3@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_37s45m16h_06-06-19_avg1_1D_raw.SGL' 0
%     '190606/resolution27umPlanar_BK31[CNT]_4@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_23s48m16h_06-06-19_avg1_1D_raw.SGL' 0
%     '190606/resolution27umPlanar_BK31[CNT]_5@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_16s50m16h_06-06-19_avg1_1D_raw.SGL' 0
%     '190606/resolution27umPlanar_BK31[CNT]_6@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_03s52m16h_06-06-19_avg1_1D_raw.SGL' 0
%     '190606/resolution27umPlanar_BK31[CNT]_7@0nm_t0[-1000]_dx[0µm]_dy[20µm]_dt[4ns]_49s53m16h_06-06-19_avg1_1D_raw.SGL' 4e-6
%     '190606/resolution27umPlanar_BK31[CNT]_8@0nm_t0[-1000]_dx[0µm]_dy[20µm]_dt[4ns]_47s55m16h_06-06-19_avg1_1D_raw.SGL' 4e-6
%     '190606/resolution27umPlanar_BK31[CNT]_9@0nm_t0[-1500]_dx[0µm]_dy[20µm]_dt[4ns]_04s58m16h_06-06-19_avg1_1D_raw.SGL' 6e-6
%     '190606/resolution27umPlanar_BK31[CNT]_10@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_57s59m16h_06-06-19_avg1_1D_raw.SGL' 7e-6
%     '190606/resolution27umPlanar_BK31[CNT]_11@0nm_t0[-2000]_dx[0µm]_dy[20µm]_dt[4ns]_59s01m17h_06-06-19_avg1_1D_raw.SGL' 8e-6
%     '190606/resolution27umPlanar_BK31[CNT]_12@0nm_t0[-2500]_dx[0µm]_dy[20µm]_dt[4ns]_43s03m17h_06-06-19_avg1_1D_raw.SGL' 10e-6
%     '190606/resolution27umPlanar_BK31[CNT]_13@0nm_t0[-2750]_dx[0µm]_dy[20µm]_dt[4ns]_56s05m17h_06-06-19_avg1_1D_raw.SGL' 11e-6
%     '190606/resolution27umPlanar_BK31[CNT]_14@0nm_t0[-2750]_dx[0µm]_dy[20µm]_dt[4ns]_45s07m17h_06-06-19_avg1_1D_raw.SGL' 11e-6
%     '190606/resolution27umPlanar_BK31[CNT]_15@0nm_t0[-3250]_dx[0µm]_dy[20µm]_dt[4ns]_46s09m17h_06-06-19_avg1_1D_raw.SGL' 13e-6
%     '190606/resolution27umPlanar_BK31[CNT]_16@0nm_t0[-3750]_dx[0µm]_dy[20µm]_dt[4ns]_53s11m17h_06-06-19_avg1_1D_raw.SGL' 15e-6
%     '190606/resolution27umPlanar_BK31[CNT]_17@0nm_t0[-4000]_dx[0µm]_dy[20µm]_dt[4ns]_42s13m17h_06-06-19_avg1_1D_raw.SGL' 16e-6
%     '190606/resolution27umPlanar_BK31[CNT]_18@0nm_t0[-4000]_dx[0µm]_dy[20µm]_dt[4ns]_25s15m17h_06-06-19_avg1_1D_raw.SGL' 16e-6
%     };
% 
% file_names      = diff_line_scans(:,1);
% trigger_delays  = diff_line_scans(:,2);
% 
% % need to reshape to be able to loop through
% file_names     = reshape(file_names,[1 length(file_names)]);
% trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);
% 
% c0 = 1484;
% 
% samples_cut_off = 50;
% samples_t0_correct = -12;


%% 190613 resolution27umPlanar2 - BK31[CNT]

% dim = 3;
% 
% file_name = '190613/resolution27umPlanar2_BK31[CNT]_2D@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_03s11m17h_13-06-19_avg1_2D_raw.SGL';
% 
% trigger_delay = 0;

% dim = 2;
% 
% diff_line_scans = {
%     '190613/resolution27umPlanar2_BK31[CNT]_1@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_00s28m17h_13-06-19_avg1_1D_raw.SGL' 0
%     '190613/resolution27umPlanar2_BK31[CNT]_2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_00s30m17h_13-06-19_avg1_1D_raw.SGL' 0
%     '190613/resolution27umPlanar2_BK31[CNT]_3@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_41s31m17h_13-06-19_avg1_1D_raw.SGL' 0
%     '190613/resolution27umPlanar2_BK31[CNT]_4@0nm_t0[-250]_dx[0µm]_dy[20µm]_dt[4ns]_42s33m17h_13-06-19_avg1_1D_raw.SGL' 1e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_5@0nm_t0[-250]_dx[0µm]_dy[20µm]_dt[4ns]_53s35m17h_13-06-19_avg1_1D_raw.SGL' 1e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_6@0nm_t0[-500]_dx[0µm]_dy[20µm]_dt[4ns]_47s37m17h_13-06-19_avg1_1D_raw.SGL' 2e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_7@0nm_t0[-750]_dx[0µm]_dy[20µm]_dt[4ns]_54s39m17h_13-06-19_avg1_1D_raw.SGL' 3e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_8@0nm_t0[-1000]_dx[0µm]_dy[20µm]_dt[4ns]_09s42m17h_13-06-19_avg1_1D_raw.SGL' 4e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_9@0nm_t0[-1250]_dx[0µm]_dy[20µm]_dt[4ns]_17s44m17h_13-06-19_avg1_1D_raw.SGL' 5e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_10@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_56s46m17h_13-06-19_avg1_1D_raw.SGL' 7e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_11@0nm_t0[-2000]_dx[0µm]_dy[20µm]_dt[4ns]_22s49m17h_13-06-19_avg1_1D_raw.SGL' 8e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_12@0nm_t0[-2500]_dx[0µm]_dy[20µm]_dt[4ns]_30s51m17h_13-06-19_avg1_1D_raw.SGL' 10e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_13@0nm_t0[-2500]_dx[0µm]_dy[20µm]_dt[4ns]_43s53m17h_13-06-19_avg1_1D_raw.SGL' 10e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_14@0nm_t0[-2750]_dx[0µm]_dy[20µm]_dt[4ns]_00s56m17h_13-06-19_avg1_1D_raw.SGL' 11e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_15@0nm_t0[-3250]_dx[0µm]_dy[20µm]_dt[4ns]_10s58m17h_13-06-19_avg1_1D_raw.SGL' 13e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_16@0nm_t0[-3500]_dx[0µm]_dy[20µm]_dt[4ns]_22s00m18h_13-06-19_avg1_1D_raw.SGL' 14e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_17@0nm_t0[-4000]_dx[0µm]_dy[20µm]_dt[4ns]_21s02m18h_13-06-19_avg1_1D_raw.SGL' 16e-6
%     '190613/resolution27umPlanar2_BK31[CNT]_18@0nm_t0[-4000]_dx[0µm]_dy[20µm]_dt[4ns]_27s04m18h_13-06-19_avg1_1D_raw.SGL' 16e-6
%     };
% 
% file_names      = diff_line_scans(:,1);
% trigger_delays  = diff_line_scans(:,2);
% 
% % need to reshape to be able to loop through
% file_names     = reshape(file_names,[1 length(file_names)]);
% trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);
% 
% c0 = 1484;
% 
% samples_cut_off = 50;
% samples_t0_correct = -12;


%% 190618 resolution27umPlanar3 (flipped in xy) - BK31[CNT]

% dim = 3;
% 
% file_name = '190618/resolution27umPlanar3_BK31[CNT]_2D@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_44s26m12h_18-06-19_avg1_2D_raw.SGL';
% 
% trigger_delay = 0;

% dim = 2;
% 
% diff_line_scans = {
%     '190618/resolution27umPlanar3_BK31[CNT]_1@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_39s34m12h_18-06-19_avg1_1D_raw.SGL' 0
%     '190618/resolution27umPlanar3_BK31[CNT]_2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_34s36m12h_18-06-19_avg1_1D_raw.SGL' 0
%     '190618/resolution27umPlanar3_BK31[CNT]_3@0nm_t0[-250]_dx[0µm]_dy[20µm]_dt[4ns]_24s38m12h_18-06-19_avg1_1D_raw.SGL' 1e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_4@0nm_t0[-250]_dx[0µm]_dy[20µm]_dt[4ns]_05s40m12h_18-06-19_avg1_1D_raw.SGL' 1e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_5@0nm_t0[-500]_dx[0µm]_dy[20µm]_dt[4ns]_43s41m12h_18-06-19_avg1_1D_raw.SGL' 2e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_6@0nm_t0[-750]_dx[0µm]_dy[20µm]_dt[4ns]_22s43m12h_18-06-19_avg1_1D_raw.SGL' 3e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_7@0nm_t0[-750]_dx[0µm]_dy[20µm]_dt[4ns]_10s45m12h_18-06-19_avg1_1D_raw.SGL' 3e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_8@0nm_t0[-1000]_dx[0µm]_dy[20µm]_dt[4ns]_44s46m12h_18-06-19_avg1_1D_raw.SGL' 4e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_9@0nm_t0[-1250]_dx[0µm]_dy[20µm]_dt[4ns]_02s49m12h_18-06-19_avg1_1D_raw.SGL' 5e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_10@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_50s50m12h_18-06-19_avg1_1D_raw.SGL' 7e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_11@0nm_t0[-2000]_dx[0µm]_dy[20µm]_dt[4ns]_29s52m12h_18-06-19_avg1_1D_raw.SGL' 8e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_12@0nm_t0[-2500]_dx[0µm]_dy[20µm]_dt[4ns]_13s54m12h_18-06-19_avg1_1D_raw.SGL' 10e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_13@0nm_t0[-2500]_dx[0µm]_dy[20µm]_dt[4ns]_51s55m12h_18-06-19_avg1_1D_raw.SGL' 10e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_14@0nm_t0[-2750]_dx[0µm]_dy[20µm]_dt[4ns]_38s57m12h_18-06-19_avg1_1D_raw.SGL' 11e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_15@0nm_t0[-3250]_dx[0µm]_dy[20µm]_dt[4ns]_29s59m12h_18-06-19_avg1_1D_raw.SGL' 13e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_16@0nm_t0[-3500]_dx[0µm]_dy[20µm]_dt[4ns]_13s01m13h_18-06-19_avg1_1D_raw.SGL' 14e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_17@0nm_t0[-4000]_dx[0µm]_dy[20µm]_dt[4ns]_22s03m13h_18-06-19_avg1_1D_raw.SGL' 16e-6
%     '190618/resolution27umPlanar3_BK31[CNT]_18@0nm_t0[-4000]_dx[0µm]_dy[20µm]_dt[4ns]_06s05m13h_18-06-19_avg1_1D_raw.SGL' 16e-6
%     };
% 
% file_names      = diff_line_scans(:,1);
% trigger_delays  = diff_line_scans(:,2);
% 
% % need to reshape to be able to loop through
% file_names     = reshape(file_names,[1 length(file_names)]);
% trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);
% 
% c0 = 1484;
% 
% samples_cut_off = 50;
% samples_t0_correct = -12;


%% 190618 resolution27umPlanar4 (flipped sensor xy) - BK31[CNT]

% dim = 3;
% 
% file_name = '190618/resolution27umPlanar4_BK31[CNT]_2D@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_03s49m17h_18-06-19_avg1_2D_raw.SGL';
% 
% trigger_delay = 0;

% dim = 2;
% 
% diff_line_scans = {
%     '190618/resolution27umPlanar4_BK31[CNT]_1@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_13s56m17h_18-06-19_avg1_1D_raw.SGL' 0
%     '190618/resolution27umPlanar4_BK31[CNT]_2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_58s57m17h_18-06-19_avg1_1D_raw.SGL' 0
%     '190618/resolution27umPlanar4_BK31[CNT]_3@0nm_t0[-250]_dx[0µm]_dy[20µm]_dt[4ns]_54s59m17h_18-06-19_avg1_1D_raw.SGL' 1e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_4@0nm_t0[-250]_dx[0µm]_dy[20µm]_dt[4ns]_35s01m18h_18-06-19_avg1_1D_raw.SGL' 1e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_5@0nm_t0[-500]_dx[0µm]_dy[20µm]_dt[4ns]_33s03m18h_18-06-19_avg1_1D_raw.SGL' 2e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_6@0nm_t0[-750]_dx[0µm]_dy[20µm]_dt[4ns]_31s05m18h_18-06-19_avg1_1D_raw.SGL' 3e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_7@0nm_t0[-750]_dx[0µm]_dy[20µm]_dt[4ns]_31s07m18h_18-06-19_avg1_1D_raw.SGL' 3e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_8@0nm_t0[-1000]_dx[0µm]_dy[20µm]_dt[4ns]_10s09m18h_18-06-19_avg1_1D_raw.SGL' 4e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_9@0nm_t0[-1250]_dx[0µm]_dy[20µm]_dt[4ns]_18s11m18h_18-06-19_avg1_1D_raw.SGL' 5e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_10@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_37s13m18h_18-06-19_avg1_1D_raw.SGL' 7e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_11@0nm_t0[-2000]_dx[0µm]_dy[20µm]_dt[4ns]_17s15m18h_18-06-19_avg1_1D_raw.SGL' 8e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_12@0nm_t0[-2500]_dx[0µm]_dy[20µm]_dt[4ns]_09s17m18h_18-06-19_avg1_1D_raw.SGL' 10e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_13@0nm_t0[-2500]_dx[0µm]_dy[20µm]_dt[4ns]_52s18m18h_18-06-19_avg1_1D_raw.SGL' 10e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_14@0nm_t0[-2750]_dx[0µm]_dy[20µm]_dt[4ns]_08s21m18h_18-06-19_avg1_1D_raw.SGL' 11e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_15@0nm_t0[-3250]_dx[0µm]_dy[20µm]_dt[4ns]_53s22m18h_18-06-19_avg1_1D_raw.SGL' 13e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_16@0nm_t0[-3500]_dx[0µm]_dy[20µm]_dt[4ns]_42s24m18h_18-06-19_avg1_1D_raw.SGL' 14e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_17@0nm_t0[-4000]_dx[0µm]_dy[20µm]_dt[4ns]_47s27m18h_18-06-19_avg1_1D_raw.SGL' 16e-6
%     '190618/resolution27umPlanar4_BK31[CNT]_18@0nm_t0[-4000]_dx[0µm]_dy[20µm]_dt[4ns]_47s29m18h_18-06-19_avg1_1D_raw.SGL' 16e-6
%     };
% 
% file_names      = diff_line_scans(:,1);
% trigger_delays  = diff_line_scans(:,2);
% 
% % need to reshape to be able to loop through
% file_names     = reshape(file_names,[1 length(file_names)]);
% trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);
% 
% c0 = 1484;
% 
% samples_cut_off = 50;
% samples_t0_correct = -12;


%% 190709 resolution27umPlanar with tip-tilted excitation beam at 5 angles

% dim = 2;
% 
% diff_line_scans = {
%     '190709/resolution27umPlanar_BK31[CNT]_10_rad1@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_25s19m21h_09-07-19_avg1_1D_raw.SGL' 7e-6
%     '190709/resolution27umPlanar_BK31[CNT]_10_rad2@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_06s25m21h_09-07-19_avg1_1D_raw.SGL' 7e-6
%     '190709/resolution27umPlanar_BK31[CNT]_10_rad3@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_08s29m21h_09-07-19_avg1_1D_raw.SGL' 7e-6
%     '190709/resolution27umPlanar_BK31[CNT]_10_rad4@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_56s31m21h_09-07-19_avg1_1D_raw.SGL' 7e-6
%     '190709/resolution27umPlanar_BK31[CNT]_10_rad5@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_40s35m21h_09-07-19_avg1_1D_raw.SGL' 7e-6
%     };
% 
% file_names      = diff_line_scans(:,1);
% trigger_delays  = diff_line_scans(:,2);
% 
% % need to reshape to be able to loop through
% file_names     = reshape(file_names,[1 length(file_names)]);
% trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);
% 
% c0 = 1484;
% 
% samples_cut_off = 50;
% samples_t0_correct = -12;


%% 190830 cause for resolution asymmetry (phantom at 6 mm with central/off-centre excitation beam)

% file_name = '190830\resolution27umPlanar_ND1_BK31[CNT]_1D_central@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_11s05m16h_30-08-19_avg1_1D_raw.SGL';
% % file_name = '190830\resolution27umPlanar_ND1_BK31[CNT]_1D_offcentre1@0nm_t0[-1750]_dx[0µm]_dy[20µm]_dt[4ns]_58s58m16h_30-08-19_avg1_1D_raw.SGL';
% 
% c0 = 1484;
% 
% trigger_delay = 7e-6;
% 
% dim = 2;
% 
% samples_cut_off = 50;
% samples_t0_correct = -9;


%% 190927 resolution27umPlanar up to 10 mm on trolley scanner with straightish fibre

% dim = 3;
% 
% file_name = '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos0@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_22s15m14h_27-09-19_avg1_2D_raw.SGL';
% 
% trigger_delay = 0;

% dim = 2;
% 
% diff_line_scans = {
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos0@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_17s05m18h_27-09-19_avg1_1D_raw.SGL' 0
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos1@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_16s09m18h_27-09-19_avg1_1D_raw.SGL' 0
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_05s11m18h_27-09-19_avg1_1D_raw.SGL' 0
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos3@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_48s13m18h_27-09-19_avg1_1D_raw.SGL' 1e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos4@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_28s17m18h_27-09-19_avg1_1D_raw.SGL' 1e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos5@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_03s21m18h_27-09-19_avg1_1D_raw.SGL' 2e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos6@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_12s23m18h_27-09-19_avg1_1D_raw.SGL' 2e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos7@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_45s25m18h_27-09-19_avg1_1D_raw.SGL' 3e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos8@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_51s28m18h_27-09-19_avg1_1D_raw.SGL' 4e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos9@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_52s31m18h_27-09-19_avg1_1D_raw.SGL' 4e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos10@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_58s34m18h_27-09-19_avg1_1D_raw.SGL' 5e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos11@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_20s38m18h_27-09-19_avg1_1D_raw.SGL' 7e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos12@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_07s41m18h_27-09-19_avg1_1D_raw.SGL' 8e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos13@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_52s43m18h_27-09-19_avg1_1D_raw.SGL' 9e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos14@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_42s46m18h_27-09-19_avg1_1D_raw.SGL' 10e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos15@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_09s50m18h_27-09-19_avg1_1D_raw.SGL' 12e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos16@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_29s52m18h_27-09-19_avg1_1D_raw.SGL' 13e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos17@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_19s55m18h_27-09-19_avg1_1D_raw.SGL' 15e-6
%     '190927/resolution27umPlanar_BK31[CNT]_trolley_straight2_2D_pos18@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_43s58m18h_27-09-19_avg1_1D_raw.SGL' 16e-6
%     };
% 
% file_names      = diff_line_scans(:,1);
% trigger_delays  = diff_line_scans(:,2);
% 
% % need to reshape to be able to loop through
% file_names     = reshape(file_names,[1 length(file_names)]);
% trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);
% 
% c0 = 1484;
% 
% samples_cut_off = 50;
% samples_t0_correct = -14; % starts at 13


%% 191029 resolution27umPlanar up to 10 mm on trolley scanner with scrambled fibre

% dim = 3;
% 
% file_name = '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_2D@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_07s23m14h_29-10-19_avg1_2D_raw.SGL';
% 
% trigger_delay = 0;

% dim = 2;
% 
% diff_line_scans = {
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos1@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_41s32m14h_29-10-19_avg1_1D_raw.SGL' 0
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_04s35m14h_29-10-19_avg1_1D_raw.SGL' 0
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos3@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_05s37m14h_29-10-19_avg1_1D_raw.SGL' 0
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos4@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_49s39m14h_29-10-19_avg1_1D_raw.SGL' 1e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos5@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_31s41m14h_29-10-19_avg1_1D_raw.SGL' 2e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos6@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_21s44m14h_29-10-19_avg1_1D_raw.SGL' 3e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos7@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_10s46m14h_29-10-19_avg1_1D_raw.SGL' 4e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos8@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_03s48m14h_29-10-19_avg1_1D_raw.SGL' 4e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos9@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_12s50m14h_29-10-19_avg1_1D_raw.SGL' 5e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos10@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_14s52m14h_29-10-19_avg1_1D_raw.SGL' 6e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos11@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_22s54m14h_29-10-19_avg1_1D_raw.SGL' 8e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos12@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_19s56m14h_29-10-19_avg1_1D_raw.SGL' 9e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos13@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_04s58m14h_29-10-19_avg1_1D_raw.SGL' 10e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos14@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_01s00m15h_29-10-19_avg1_1D_raw.SGL' 12e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos15@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_44s02m15h_29-10-19_avg1_1D_raw.SGL' 13e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos16@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_41s04m15h_29-10-19_avg1_1D_raw.SGL' 14e-6
%     '191029/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos17@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_36s06m15h_29-10-19_avg1_1D_raw.SGL' 16e-6
%     };
% 
% file_names      = diff_line_scans(:,1);
% trigger_delays  = diff_line_scans(:,2);
% 
% % need to reshape to be able to loop through
% file_names     = reshape(file_names,[1 length(file_names)]);
% trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);
% 
% c0 = 1484;
% 
% samples_cut_off = 50;
% samples_t0_correct = -14;


%% 191031 resolution27umPlanar up to 10 mm on trolley scanner with scrambled fibre & (more) centralised phantom

% dim = 3;
% 
% file_name = '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_2D@0nm_t0[0]_dx[100µm]_dy[100µm]_dt[4ns]_27s34m15h_31-10-19_avg1_2D_raw.SGL';
% 
% trigger_delay = 0;

dim = 2;

diff_line_scans = {
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos1@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_37s42m15h_31-10-19_avg1_1D_raw.SGL' 0
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos2@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_29s44m15h_31-10-19_avg1_1D_raw.SGL' 0
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos3@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_38s46m15h_31-10-19_avg1_1D_raw.SGL' 1e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos4@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_15s48m15h_31-10-19_avg1_1D_raw.SGL' 1e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos5@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_48s50m15h_31-10-19_avg1_1D_raw.SGL' 2e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos6@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_37s52m15h_31-10-19_avg1_1D_raw.SGL' 3e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos7@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_36s54m15h_31-10-19_avg1_1D_raw.SGL' 3e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos8@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_26s56m15h_31-10-19_avg1_1D_raw.SGL' 4e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos9@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_13s58m15h_31-10-19_avg1_1D_raw.SGL' 5e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos10@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_00s00m16h_31-10-19_avg1_1D_raw.SGL' 6e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos11@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_43s02m16h_31-10-19_avg1_1D_raw.SGL' 8e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos12@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_40s04m16h_31-10-19_avg1_1D_raw.SGL' 9e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos13@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_31s06m16h_31-10-19_avg1_1D_raw.SGL' 10e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos14@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_55s08m16h_31-10-19_avg1_1D_raw.SGL' 12e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos15@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_38s10m16h_31-10-19_avg1_1D_raw.SGL' 13e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos16@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_34s12m16h_31-10-19_avg1_1D_raw.SGL' 14e-6
    '191031/resolution27umPlanar_BK31[CNT]_trolley_scrambled_1D_pos17@0nm_t0[0]_dx[0µm]_dy[20µm]_dt[4ns]_36s14m16h_31-10-19_avg1_1D_raw.SGL' 16e-6
    };

file_names      = diff_line_scans(:,1);
trigger_delays  = diff_line_scans(:,2);

% need to reshape to be able to loop through
file_names     = reshape(file_names,[1 length(file_names)]);
trigger_delays = reshape(trigger_delays,[1 length(trigger_delays)]);

c0 = 1484;

samples_cut_off = 50;
samples_t0_correct = -14;

