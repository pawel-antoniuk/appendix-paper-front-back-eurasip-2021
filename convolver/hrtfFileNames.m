% Sławomir Zieliński 2021
% Bialystok University of Technology

function [hrtfName] = hrtfFileNames()
% SOFA HRTF Filenames
% Outputs a stucture
hrtfName = {};
hrtfName{1} = 'mit\mit_kemar_normal_pinna_48kHz_sz.sofa';
hrtfName{2} = 'mit\mit_kemar_large_pinna_48kHz_sz.sofa';
hrtfName{3} = 'sadie\D1_48K_24bit_256tap_FIR_SOFA.sofa';
hrtfName{4} = 'sadie\D2_48K_24bit_256tap_FIR_SOFA.sofa';
hrtfName{5} = 'sadie\H3_48K_24bit_256tap_FIR_SOFA.sofa';
hrtfName{6} = 'sadie\H4_48K_24bit_256tap_FIR_SOFA.sofa';
hrtfName{7} = 'sadie\H5_48K_24bit_256tap_FIR_SOFA.sofa';
hrtfName{8} = 'sadie\H6_48K_24bit_256tap_FIR_SOFA.sofa';
hrtfName{9} = 'sadie\H7_48K_24bit_256tap_FIR_SOFA.sofa';
hrtfName{10} = 'listen\LISTEN_1003_IRC_1003_C_HRIR.sofa';
hrtfName{11} = 'listen\LISTEN_1046_IRC_1046_C_HRIR.sofa';
hrtfName{12} = 'listen\LISTEN_1054_IRC_1054_C_HRIR.sofa';
hrtfName{13} = 'listen\LISTEN_1006_IRC_1006_C_HRIR.sofa';
hrtfName{14} = 'listen\LISTEN_1002_IRC_1002_C_HRIR.sofa';
hrtfName{15} = 'listen\LISTEN_1004_IRC_1004_C_HRIR.sofa';
hrtfName{16} = 'listen\LISTEN_1005_IRC_1005_C_HRIR.sofa';
hrtfName{17} = 'cipic\subject_012_48kHz_sz.sofa';
hrtfName{18} = 'cipic\subject_015_48kHz_sz.sofa';
hrtfName{19} = 'cipic\subject_020_48kHz_sz.sofa';
hrtfName{20} = 'cipic\subject_028_48kHz_sz.sofa';
hrtfName{21} = 'cipic\subject_051_48kHz_sz.sofa';
hrtfName{22} = 'cipic\subject_147_48kHz_sz.sofa';
hrtfName{23} = 'cipic\subject_148_48kHz_sz.sofa';
hrtfName{24} = 'th-koln\HRIR_L2702_NF050.sofa';
hrtfName{25} = 'th-koln\HRIR_L2702_NF075.sofa';
hrtfName{26} = 'th-koln\HRIR_L2702_NF100.sofa';
hrtfName{27} = 'th-koln\HRIR_L2702_NF150.sofa';
hrtfName{28} = 'th-koln\HMSII.sofa';
hrtfName{29} = 'th-koln\HMSII_Cap.sofa';
hrtfName{30} = 'th-koln\HMSII_OculusRift.sofa';
hrtfName{31} = 'tu-berlin\FABIAN_HRIR_measured_HATO_0_48kHz_sz.sofa';
hrtfName{32} = 'tu-berlin\qu_kemar_anechoic_0.5m_48kHz_sz.sofa';
hrtfName{33} = 'tu-berlin\qu_kemar_anechoic_1m_48kHz_sz.sofa';
hrtfName{34} = 'tu-berlin\qu_kemar_anechoic_2m_48kHz_sz.sofa';
hrtfName{35} = 'tu-berlin\qu_kemar_anechoic_3m_48kHz_sz.sofa';
hrtfName{36} = 'riec\RIEC_hrir_subject_001.sofa';
hrtfName{37} = 'riec\RIEC_hrir_subject_002.sofa';
hrtfName{38} = 'riec\RIEC_hrir_subject_003.sofa';
hrtfName{39} = 'riec\RIEC_hrir_subject_004.sofa';
hrtfName{40} = 'riec\RIEC_hrir_subject_005.sofa';
hrtfName{41} = 'riec\RIEC_hrir_subject_046.sofa'; % SAMRAI
hrtfName{42} = 'riec\RIEC_hrir_subject_080.sofa'; % KEMAR
hrtfName{43} = 'ari\hrtf b_nh2.sofa';
hrtfName{44} = 'ari\hrtf b_nh4.sofa';
hrtfName{45} = 'ari\hrtf b_nh5.sofa';
hrtfName{46} = 'ari\hrtf b_nh8.sofa';
hrtfName{47} = 'ari\hrtf b_nh10.sofa';
hrtfName{48} = 'ari\hrtf b_nh169.sofa';
hrtfName{49} = 'ari\hrtf b_nh172.sofa';
hrtfName{50} = 'scut\SCUT_KEMAR_radius_1_48kHz_sz.sofa';
hrtfName{51} = 'scut\SCUT_KEMAR_radius_0.5_48kHz_sz.sofa';
hrtfName{52} = 'clubfritz\ClubFritz1_48kHz_sz.sofa';
hrtfName{53} = 'clubfritz\ClubFritz3_48kHz_sz.sofa';
hrtfName{54} = 'clubfritz\ClubFritz4_48kHz_sz.sofa';
hrtfName{55} = 'clubfritz\ClubFritz6.sofa';
hrtfName{56} = 'clubfritz\ClubFritz7.sofa';
hrtfName{57} = 'clubfritz\ClubFritz8.sofa';
hrtfName{58} = 'clubfritz\ClubFritz9.sofa';
hrtfName{59} = 'viking\subj_A.sofa';
hrtfName{60} = 'viking\subj_B.sofa';
hrtfName{61} = 'viking\subj_C.sofa';
hrtfName{62} = 'viking\subj_D.sofa';
hrtfName{63} = 'viking\subj_E.sofa';
hrtfName{64} = 'viking\subj_F.sofa';
hrtfName{65} = 'viking\subj_G.sofa';
hrtfName{66} = 'aachen\HRTF_perFrequencyInterpolation.sofa';
hrtfName{67} = 'aachen\Kemar_HRTF_sofa.sofa';
hrtfName{68} = 'hutubs\pp1_HRIRs_measured_48kHz_sz.sofa';
hrtfName{69} = 'hutubs\pp2_HRIRs_measured_48kHz_sz.sofa';
hrtfName{70} = 'hutubs\pp3_HRIRs_measured_48kHz_sz.sofa';
hrtfName{71} = 'hutubs\pp4_HRIRs_measured_48kHz_sz.sofa';
hrtfName{72} = 'hutubs\pp5_HRIRs_measured_48kHz_sz.sofa';
hrtfName{73} = 'hutubs\pp6_HRIRs_measured_48kHz_sz.sofa';
hrtfName{74} = 'hutubs\pp7_HRIRs_measured_48kHz_sz.sofa';
end















