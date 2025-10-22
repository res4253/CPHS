function dx = linear_dynamics(~,x,u,p)
            %%
            G_L = x(1);
            X_L = x(2);
            Q_L1 = x(3);%%
            Q_L2 = x(4);%%
            I_L1 = x(5);
            I_L2 = x(6);
            I_L = x(7);  

            %%
            u_m = u(1);
            u_i = u(2);
            mode = u(3);%%
            
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
        
            R_a = f_c*k_abs*Q_L2;

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

            %%
            xe = [G_L-Ge; X_L-Xe; I_L1-i_b/k_d; I_L2-i_b/k_d; I_L-i_b/k_cl];
            dxe = A*xe + B_1*(u_i-i_b) + B_2*R_a;

            dQ_L1 = -k_tau*Q_L1 + u_m;
            dQ_L2 = -k_abs*Q_L2 + k_tau*Q_L1;

            %%
            dx = [dxe(1); dxe(2); dQ_L1; dQ_L2; dxe(3); dxe(4); dxe(5)];
        end