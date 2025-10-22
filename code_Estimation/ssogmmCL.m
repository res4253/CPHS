classdef ssogmmCL
    properties 
        n_x
        n_p
        n_u
        C
        p_0 % パラメータの公称値
        x_0
        u_0
        u_inb
        BW
        S_g
        V_g
        V_I
        p_2
        f_c
        G_b
        k_tau
        k_abs_0
        k_abs_1
        k_d
        k_cl
        I_b
        TDI_basal
        TDI_whole
        S_I
        ubX
        ubP
        lbX
        lbP
        p_c
    end
    methods
        function self = ssogmmCL(options)
            arguments
                options.u_inb (1,1) {mustBeNumeric} = NaN;
                options.BW (1,1) {mustBeNumeric} = NaN;
                options.S_g (1,1) {mustBeNumeric} = NaN; 
                options.V_g (1,1) {mustBeNumeric} = NaN; 
                options.V_I (1,1) {mustBeNumeric} = NaN; 
                options.p_2 (1,1) {mustBeNumeric} = NaN; 
                options.f_c (1,1) {mustBeNumeric} = NaN;   
                options.G_b (1,1) {mustBeNumeric} = NaN;
                options.k_tau (1,1) {mustBeNumeric} = NaN;
                options.k_abs_0 (1,1) {mustBeNumeric} = NaN;
                options.k_abs_1 (1,1) {mustBeNumeric} = NaN;
                options.k_d (1,1) {mustBeNumeric} = NaN;
                options.k_cl (1,1) {mustBeNumeric} = NaN; 
                options.ubX (:,1) {mustBeNumeric} = NaN;
                options.lbX (:,1) {mustBeNumeric} = NaN;
                options.ubP (:,1) {mustBeNumeric} = NaN;
                options.lbP (:,1) {mustBeNumeric} = NaN;
                options.p_c = [];
            end          

            load("t1dms_data/nominal_SOGMM.mat", "nominal_SOGMM");

            self.n_x = 7;
            self.n_u = 3;
            self.n_p = 14;
        
            self.u_inb = nominal_SOGMM.u_inb;
            self.BW = nominal_SOGMM.BW;
            self.S_g = nominal_SOGMM.S_g; 
            self.V_g = nominal_SOGMM.V_g; 
            self.V_I = nominal_SOGMM.V_I; 
            self.p_2 = nominal_SOGMM.p_2; 
            self.f_c = nominal_SOGMM.f_c; 
            self.G_b = nominal_SOGMM.G_b; 
            self.k_tau = nominal_SOGMM.k_tau;
            self.k_abs_0 = nominal_SOGMM.k_abs*0.8;
            self.k_abs_1 = nominal_SOGMM.k_abs*1.3;
            self.k_d = nominal_SOGMM.k_d; 
            self.k_cl = nominal_SOGMM.k_cl; 
            self.ubX = Inf(self.n_x,1);
            self.lbX = zeros(self.n_x,1);
            self.ubP = Inf(self.n_p,1);
            self.lbP = zeros(self.n_p,1);
            self.p_c = options.p_c;

            if ~isempty(self.p_c) 
                if slength(self.p_c) ~= self.n_p           
                    error("p_c の次元が p の次元と一致しません")
                end
            else
                self.p_c = nan(self.n_p,1);
            end

            % option で入力された場合は selfに代入する 
            fields = fieldnames(options);
            for i = 1:numel(fields)
                field = fields{i};
                if ~(isnan(options.(field)))
                    self.(field) = options.(field);
                end
            end        

            self.I_b = self.u_inb/(self.k_cl*self.V_I*self.BW); 
            self.TDI_basal = (self.u_inb*1440)/1000; % U scale
            self.TDI_whole = 38.9;
            self.S_I = exp(- 6.4417 - 0.063546*self.TDI_whole + 0.057944*self.TDI_basal);  

            self.p_0 = [
                self.G_b;
                self.V_I;
                self.S_I;
                self.k_tau;
                self.k_abs_0;
                self.k_abs_1;
                self.k_d;
                self.k_cl;
                self.S_g;
                self.V_g;
                self.p_2;
                self.BW;
                self.f_c;
                self.I_b
                ];
            self.x_0 = [
                self.G_b; 
                0; 
                0; 
                0; 
                self.u_inb/self.k_d; 
                self.u_inb/self.k_d; 
                self.u_inb/self.k_cl
                ];            
            self.u_0 = [
                0
                self.u_inb
                ];

            self.C = [1 0 0 0 0 0 0];            

        end
        function dx = dynamics(~,~,x,u,p)
            G = x(1);
            X_I = x(2);
            Q_1 = x(3);
            Q_2 = x(4);
            I_sc1 = x(5);
            I_sc2 = x(6);
            I_p = x(7);
        
            u_m = u(1);
            u_i = u(2);
            mode = u(3);
        
            if mode == 0
                k_abs = p(5);
            else
                k_abs = p(6);
            end
        
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
        
            dG = -(S_g + X_I)*G + S_g*G_b + R_a/(BW*V_g);
            dX_I = -p_2*X_I + p_2*S_I*(I_p/(BW*V_I) - I_b);
            dQ_1 = -k_tau*Q_1 + u_m;
            dQ_2 = -k_abs*Q_2 + k_tau*Q_1;
            dI_sc1 = -k_d*I_sc1 + u_i;
            dI_sc2 = -k_d*I_sc2 + k_d*I_sc1;
            dI_p = -k_cl*I_p + k_d*I_sc2;
        
            dx = [dG; dX_I; dQ_1; dQ_2; dI_sc1; dI_sc2; dI_p];
        end
        function dx = linear_dynamics(~,~,x,u,p)

            

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
    end
end