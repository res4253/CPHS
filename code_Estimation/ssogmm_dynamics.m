function dx = ssogmm_dynamics(~,x,u,p)
            %%
            G = x(1);
            X_I = x(2);
            Q_1 = x(3);
            Q_2 = x(4);
            I_sc1 = x(5);
            I_sc2 = x(6);
            I_p = x(7);
        
            %%
            u_m = u(1);
            u_i = u(2);
            mode = u(3);
        
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
        
            R_a = f_c*k_abs*Q_2;
        
            %%
            dG = -(S_g + X_I)*G + S_g*G_b + R_a/(BW*V_g);
            dX_I = -p_2*X_I + p_2*S_I*(I_p/(BW*V_I) - I_b);
            dQ_1 = -k_tau*Q_1 + u_m;
            dQ_2 = -k_abs*Q_2 + k_tau*Q_1;
            dI_sc1 = -k_d*I_sc1 + u_i;
            dI_sc2 = -k_d*I_sc2 + k_d*I_sc1;
            dI_p = -k_cl*I_p + k_d*I_sc2;
        
            %%
            dx = [dG; dX_I; dQ_1; dQ_2; dI_sc1; dI_sc2; dI_p];
        end