close all
% clear
clc

ratio = [];
for subject_num = 1:10

        % load(['id',num2str(subject_num),'_seed',num2str(seed)],"p_esti")
        % load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed)],"p_esti","q_esti")
        load("parametor.mat")
        load("nominal_SOGMM_13p.mat")
        
        p = adult(subject_num).mean;
        nominal = nominal_SOGMM_13p;

        nominal = [nominal(2:4);nominal(6:10)];
        p = [p(2:4);p(7:11)];

        ratio = [ratio,100*(p./nominal)];
       


end
ratio = round(ratio);

Min = min(ratio,[],"all")

Max = max(ratio,[],"all")

%ssogmm paramator
% % G_b = p(1);
% % V_I = p(2);%%
% % S_I = p(3);%
% % k_tau = p(4);
% % k_abs = p(5);
% % k_abs = p(6);
% % k_d = p(7);
% % k_cl = p(8);%%
% % S_g = p(9);
% % V_g = p(10);
% % p_2 = p(11);
% % BW = p(12);
% % f_c = p(13);
% % I_b = p(14);