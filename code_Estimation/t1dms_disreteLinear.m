close all
clear
clc

Path = ['./t1dms_data/:',...
    './t1dms_data/subjects/:',...
    './t1dms_data/simu_data/:'];
addpath(Path);


load('temp.mat') %サンプリングデータのロード


load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed),'.mat'])
p = p_esti;

%sampling time
N = size(data_ssogmm.ts,2);
N_step = 1:N;

%for obs
G = data.G.signals.values';
R = G(1:Ts:end);

Xs = data_ssogmm.xs;
%feedback
Y = G(data_ssogmm.ts);
Y = data_ssogmm.Gs;

%%input
u_m = data_ssogmm.meal;
u_i = data_ssogmm.insulin;
modes = data_ssogmm.modes;

[xe,A_d,B_d1,B_d2,C,i_b,Ge] = linear_matrix(p,Ts);

%%Ra
R_a = get_Ra(data_ssogmm.xs(4,:),p,modes,N);

%%%linear
delta_x=zeros(5,N);
x=zeros(5,N);

delta_x(:,1) = data_ssogmm.xs([1 2 5 6 7],1)-xe;
x(:,1) = data_ssogmm.xs([1 2 5 6 7],1);

for k=2:N

    delta_u_i = u_i(:,k-1)-i_b;

    delta_x(:,k) = A_d*delta_x(:,k-1) + B_d1*delta_u_i + B_d2*R_a(:,k-1);
    x(:,k) = delta_x(:,k) + xe;
end


%observer
x_ini = x_esti(:,500);
x_ini = x_ini([1 2 5 6 7]);

L = gain(p,Ts);

delta_xob = zeros(5,N);
xob = zeros(5,N);

delta_xob(:,1) = x_ini - xe;
xob(:,1) = x_ini;

for k=2:N

    delta_y = Y(:,k-1) - Ge;
    delta_u_i = u_i(:,k-1) - i_b;

    delta_xob(:,k) = A_d*delta_xob(:,k-1) + B_d1*delta_u_i + B_d2*R_a(:,k-1) + L*(delta_y - C*delta_xob(:,k-1));
    xob(:,k) = delta_xob(:,k) + xe;
end


plot(Y,'LineWidth',2)
hold on
plot(x(1,:),'.')
plot(xob(1,:),'.')
legend('uva/padova','x','xob')