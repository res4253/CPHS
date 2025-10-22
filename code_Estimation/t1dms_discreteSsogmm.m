% close all
clear
clc

Path = ['./t1dms_data/:',...
    './t1dms_data/subjects/:',...
    './t1dms_data/simu_data/:'];
addpath(Path);


%load data file
load('temp.mat') %サンプリングデータのロード

load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed),'.mat'])
p = p_esti;

%sampling time
N = size(data_ssogmm.ts,2);
N_step = 1:N;

%data
G = data.G.signals.values';
R = G(1:Ts:end);

Xs = data_ssogmm.xs;
% Xs = x_esti;
%feedback
Y = G(data_ssogmm.ts);
Y = data_ssogmm.Gs;

%%input
u_m = data_ssogmm.meal;
u_i = data_ssogmm.insulin;
modes = data_ssogmm.modes;

%Ra
R_a = get_Ra(Xs(4,:),p,modes,N);

U = [u_m;u_i;R_a;modes];


%%
Q = eps; %システム雑音
R = eps; %観測雑音
% v = randn(1,N)*sqrtm(Q);
% w = randn(1,N)*sqrtm(R);


%%ssogmm
dic = dic_ssogmmCL(p,Ts,N,U);
x = dic.dicrete(Xs,"ZOH");


%kalman
x_ini = x_esti(:,500);
x_ini = x_ini([1 2 5 6 7]);


A = dic.A;
A_di = dic.A_di;
A_dx = dic.A_dx;
C = [1 0 0 0 0];

xhat = zeros(5,N);

a = 1e-6;
P = a*eye(5);
xhat(:,1) = x_ini;

b = eye(5);
for k = 2:N

     xhatm = dic.dynamics(xhat(:,k-1),k);
     
     A_tmp = [zeros(3,2) A_di];
     A_temp =[A(xhat(1,k-1),xhat(2,k-1)) ; A_tmp];
     % Pm = A(xhat(1,k-1),xhat(2,k-1))*P*A(xhat(1,k-1),xhat(2,k-1)) + b*Q*b';
     Pm = A_temp*P*A_temp' + b*Q*b';

     g = Pm*C'/(C*Pm*C' + R);

     xhat(:,k) = xhatm + g*(Y(:,k) - C*xhatm);
     P = (eye(5) - g*C)*Pm;
end


f=figure;
plot(Y,'LineWidth',2)
hold on
plot(x(1,:),'.')
plot(xhat(1,:),'.')
legend('uva/padova','x','xhat')
