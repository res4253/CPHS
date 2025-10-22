function dx = obs(~,x,u,p,L)  
            %%
            G_ob = x(1);
            X_ob = x(2);
            Q_ob1 = x(3);%%
            Q_ob2 = x(4);%%
            I_ob1 = x(5);
            I_ob2 = x(6);
            I_ob = x(7);

            %%%
            u_m = u(1);
            u_i = u(2);
            mode = u(3);%%
            G_m = u(4);

            %%
            if mode == 0
                k_abs = p(5);
            else
                k_abs = p(6);
            end

            %%
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

            R_a = f_c*k_abs*Q_ob2;

            %%
            i_b = I_b*k_cl*V_I*BW;
            Xe = S_I*(i_b/(k_cl*BW*V_I) - I_b);
            Ge = (S_g*G_b)/(S_g+Xe);     

            %%
            A = ...
                [-(S_g+Xe), -Ge, 0, 0, 0;
                0, -p_2, 0, 0, (p_2*S_I)/(BW*V_I);
                0, 0, -k_d, 0, 0;
                0, 0, k_d, -k_d, 0;
                0, 0, 0, k_d, -k_cl
                ];
            B_1 = [0; 0; 1; 0; 0];
            B_2 = [1/(BW*V_g); 0; 0; 0; 0];
            C = [1 0 0 0 0];

            %%
            x_hat = [G_ob-Ge; X_ob-Xe; I_ob1-i_b/k_d; I_ob2-i_b/k_d; I_ob-i_b/k_cl];
            dx_hat = A*x_hat + B_1*(u_i-i_b) + B_2*R_a + L*((G_m-Ge) - C*x_hat);

            dQ_ob1 = -k_tau*Q_ob1 + u_m;
            dQ_ob2 = -k_abs*Q_ob2 + k_tau*Q_ob1;
            %%%%

            dx = [dx_hat(1:2,:); dQ_ob1 ;dQ_ob2 ; dx_hat(3:end,:)];
        end