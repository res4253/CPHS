function calculate_Ra = get_Ra(Q_2,p,modes,N)

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


R_a = zeros(1,N);
for i=1:N
    mode = modes(1,i);
    if mode == 0
        k_abs = p(5);
        R_a(1,i) = (Q_2(1,i)*k_abs*f_c)/BW;
    else
        k_abs = p(6);
        R_a(1,i) = (Q_2(1,i)*k_abs*f_c)/BW;
    end
end

calculate_Ra = R_a;
end