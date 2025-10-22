classdef sogmmCL
    properties
        n_x
        n_p
        n_u
        C
        p_0 % パラメータの公称値
        p_c
        x_0
        u_0
        y_0 
        u_inb
        BW
        S_g
        V_g
        V_I
        p_2
        f_c
        G_b
        k_tau
        k_abs
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
        J_x
        J_theta
        x_sym
        u_sym
        theta_sym
    end
    methods
        function self = sogmmCL(options) 
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
                options.k_abs (1,1) {mustBeNumeric} = NaN;
                options.k_d (1,1) {mustBeNumeric} = NaN;
                options.k_cl (1,1) {mustBeNumeric} = NaN; 
                options.ubX (:,1) {mustBeNumeric} = NaN;
                options.lbX (:,1) {mustBeNumeric} = NaN;
                options.ubP (:,1) {mustBeNumeric} = NaN;
                options.lbP (:,1) {mustBeNumeric} = NaN;
            end
            
            load("t1dms_data/nominal_SOGMM.mat", "nominal_SOGMM");

            self.n_x = 7;
            self.n_u = 2;
            self.n_p = 13; 

            self.u_inb = nominal_SOGMM.u_inb;
            self.BW = nominal_SOGMM.BW;
            self.S_g = nominal_SOGMM.S_g; 
            self.V_g = nominal_SOGMM.V_g; 
            self.V_I = nominal_SOGMM.V_I; 
            self.p_2 = nominal_SOGMM.p_2; 
            self.f_c = nominal_SOGMM.f_c; 
            self.G_b = nominal_SOGMM.G_b; 
            self.k_tau = nominal_SOGMM.k_tau;
            self.k_abs = nominal_SOGMM.k_abs;
            self.k_d = nominal_SOGMM.k_d; 
            self.k_cl = nominal_SOGMM.k_cl; 
            self.ubX = Inf(self.n_x,1);
            self.lbX = zeros(self.n_x,1);
            self.ubP = Inf(self.n_p,1);
            self.lbP = zeros(self.n_p,1);

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
            self.TDI_whole = 42;
            self.S_I = exp(- 6.4417 - 0.063546*self.TDI_whole + 0.057944*self.TDI_basal);  

            self.p_0 = [
                self.G_b;
                self.V_I;
                self.S_I;
                self.k_tau;
                self.k_abs;
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
            self.y_0 = self.G_b;

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
        
            G_b = p(1);
            V_I = p(2);
            S_I = p(3);     
            k_tau = p(4);
            k_abs = p(5);
            k_d = p(6);
            k_cl = p(7);
            S_g = p(8);
            V_g = p(9);
            p_2 = p(10);
            BW = p(11);      
            f_c = p(12);   
            I_b = p(13);            

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

        function self = get_J(self)
            self.x_sym = sym('x_sym',[1,self.n_x]);
            self.u_sym = sym('u_sym',[1,self.n_u]);
            self.theta_sym = sym('theta_sym',[1,self.n_p]);
            dx = self.dynamics(0,self.x_sym,self.u_sym,self.theta_sym);
            self.J_x = jacobian(dx,self.x_sym);
            self.J_theta = jacobian(dx,self.theta_sym);
        end

        function dxS = sensitivity_dynamics(self,t,xS,u,p)
            x = xS(1:self.n_x);
            S = reshape(xS(self.n_x+1:end),self.n_x,self.n_p);

            G = x(1);
            X_I = x(2);
            Q_1 = x(3);
            Q_2 = x(4);
            I_sc1 = x(5);
            I_sc2 = x(6);
            I_p = x(7);
        
            u_m = u(1);
            u_i = u(2);            

            G_b = p(1);
            V_I = p(2);
            S_I = p(3);     
            k_tau = p(4);
            k_abs = p(5);
            k_d = p(6);
            k_cl = p(7);
            S_g = p(8);
            V_g = p(9);
            p_2 = p(10);
            BW = p(11);      
            f_c = p(12);   
            I_b = p(13);            

            J_x = ... 
                [-S_g-X_I, -G, 0, (f_c*k_abs)/(BW*V_g), 0, 0, 0;
                0, -p_2, 0, 0, 0, 0, (S_I*p_2)/(BW*V_I);
                0, 0, -k_tau, 0, 0, 0, 0;
                0, 0, k_tau, -k_abs, 0, 0, 0;
                0, 0, 0, 0, -k_d, 0, 0;
                0, 0, 0, 0, k_d, -k_d, 0;
                0, 0, 0, 0, 0, k_d, -k_cl];


            J_theta = ...  
                [S_g, 0, 0, 0, (Q_2*f_c)/(BW*V_g), 0, 0, G_b-G, -(Q_2*f_c*k_abs)/(BW*V_g^2), 0, -(Q_2*f_c*k_abs)/(BW^2*V_g), (Q_2*k_abs)/(BW*V_g), 0; 
                0, -(I_p*S_I*p_2)/(BW*V_I^2), -p_2*(I_b - I_p/(BW*V_I)), 0, 0, 0, 0, 0, 0, -X_I-S_I*(I_b-I_p/(BW*V_I)), -(I_p*S_I*p_2)/(BW^2*V_I), 0, -S_I*p_2;
                0, 0, 0, -Q_1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, Q_1, -Q_2, 0, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, -I_sc1, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, I_sc1-I_sc2, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, I_sc2, -I_p, 0, 0, 0, 0, 0, 0];

            dx = self.dynamics(t,x,u,p);
            dS = J_x*S + J_theta;
            dxS = [dx; reshape(dS,self.n_x*self.n_p,1)];
        end
    end
end