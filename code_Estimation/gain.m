function [obserber_gain,pole] = gain(p,Ts,Q,R)

G_b = p(1);
V_I = p(2);
S_I = p(3);
k_tau = p(4);
k_d = p(7);
k_cl = p(8);
S_g = p(9);
V_g = p(10);
p_2 = p(11);
BW = p(12);
f_c = p(13);
I_b = p(14);

i_b = I_b*k_cl*V_I*BW;
Xe = S_I*(i_b/(k_cl*BW*V_I) - I_b);
Ge = (S_g*G_b)/(S_g+Xe);


A = ...
    [-(S_g+Xe), -Ge, 0, 0, 0;
    0, -p_2, 0, 0, (p_2*S_I)/(BW*V_I);
    0, 0, -k_d, 0, 0;
    0, 0, k_d, -k_d, 0;
    0, 0, 0, k_d, -k_cl
    ];
C = [1 0 0 0 0];

O = obsv(A, C);      % 可観測性行列を作成
rank_O = rank(O);     % 行列のランクを確認


if Ts ~= 0
    A = expm(A*Ts);
end


% [-0.0054;-0.0024;-0.1221;-0.0389;-0.0389]
% op = [-0.0084; -0.025; -0.123; -0.039; -0.04];
% L = place(A', C', op)';

if Ts ~= 0
    obserber_gain = dlqr(A', C', Q, R)'; % 推定ゲイン
else
    obserber_gain = lqr(A', C', Q, R)'; % 推定ゲイン
end

pole = eig(A - obserber_gain*C);

end