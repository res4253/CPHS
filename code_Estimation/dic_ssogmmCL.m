classdef dic_ssogmmCL
    properties 
        p
        Ts
        N
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
        S_I
        A
        A_di
        A_dx
        B
        B_di
        B_dx
        u_m
        u_i
        R_a
        modes
    end

    methods
        function self = dic_ssogmmCL(p,Ts,N,U)
            self.p = p;
            self.Ts = Ts;
            self.N = N;

            self.u_m = U(1,:);
            self.u_i = U(2,:);
            self.R_a = U(3,:);
            self.modes = U(4,:);

            self.G_b = p(1);
            self.V_I = p(2);
            self.S_I = p(3);
            self.k_tau = p(4);
            self.k_abs_0 = p(5);
            self.k_abs_1 = p(6);
            self.k_d = p(7);
            self.k_cl = p(8);
            self.S_g = p(9);
            self.V_g = p(10);
            self.p_2 = p(11);
            self.BW = p(12);
            self.f_c = p(13);
            self.I_b = p(14);
 

            self.A = @(Ge,Xe) ...
                [1-self.Ts*(self.S_g+Xe), -self.Ts*Ge, 0, 0, 0;
                0, 1-self.Ts*self.p_2, 0, 0, (self.Ts*self.p_2*self.S_I)/(self.BW*self.V_I)
                ];

            %%皮下インスリン動態サブモデル　行列
            A_I = [-self.k_d,0,0;
                self.k_d,-self.k_d,0;
                0,self.k_d,-self.k_cl];
            B_I = [1; 0; 0];

            self.A_di = expm(A_I*Ts);
            fun = @(tau) expm(A_I*tau);
            self.B_di = integral(fun,0,Ts,"ArrayValued",true)*B_I;

            %%状態Xについていの　行列
            A_X = -self.p_2;
            B_X = self.p_2*self.S_I;
            self.A_dx = expm(A_X*Ts);
            fun = @(tau) expm(A_X*tau);
            self.B_dx = integral(fun,0,Ts,"ArrayValued",true)*B_X;
        end
        function x = dicrete(self,xs,pattern)
            N =self.N;
            Ts = self.Ts;
            p = self.p;

            u_m = self.u_m;
            u_i = self.u_i;
            R_a = self.R_a;
            modes = self.modes;

            G = zeros(1,N);
            X = zeros(1,N);
            Q_1 = zeros(1,N);
            Q_2 = zeros(1,N);
            I_sc1 = zeros(1,N);
            I_sc2 = zeros(1,N);
            I_p = zeros(1,N);

            G(:,1) = xs(1,1);
            X(:,1) = xs(2,1);
            Q_1(:,1) = xs(3,1);
            Q_2(:,1) = xs(4,1);
            I_sc1(:,1) = xs(5,1);
            I_sc2(:,1) = xs(6,1);
            I_p(:,1) = xs(7,1);

            Q_array = [Q_1; Q_2];
            I_array = [I_sc1; I_sc2; I_p];


            if pattern == "ForwardEuler"
                for k=2:N
                    I = I_p(:,k-1)/(self.V_I*self.BW);

                    G(:,k) = G(:,k-1) + Ts*(-(self.S_g+X(:,k-1))*G(:,k-1) + self.S_g*self.G_b + R_a(:,k-1)/self.V_g);
                    X(:,k) = X(:,k-1) + Ts*(-self.p_2*X(:,k-1) + self.p_2*self.S_I*(I - self.I_b));
                    I_sc1(:,k) = I_sc1(:,k-1) + Ts*(-self.k_d*I_sc1(:,k-1) + u_i(:,k-1));
                    I_sc2(:,k) = I_sc2(:,k-1) + Ts*(-self.k_d*I_sc2(:,k-1) + self.k_d*I_sc1(:,k-1));
                    I_p(:,k) = I_p(:,k-1) + Ts*(-self.k_cl*I_p(:,k-1) + self.k_d*I_sc2(:,k-1));
                end
                x = [G; X; I_sc1; I_sc2; I_p];
                disp('ForwardEuler')
            elseif pattern == "ForwardEulerInQ"

                input = [u_m;u_i;modes];
                x = zeros(7,N);
                x(:,1) = xs(:,1);

                for k=2:N
                    x(:,k) = x(:,k-1) + Ts*ssogmm_dynamics(0,x(:,k-1),input(:,k-1),p);
                end
                disp('EulerinQ')
            elseif pattern == "ZOH"
                %ZOHによる離散化
                for k=2:N
                    I = I_array(3,k-1)/(self.V_I*self.BW);

                    G(:,k) = G(:,k-1) + Ts*(-(self.S_g+X(:,k-1))*G(:,k-1) + self.S_g*self.G_b + R_a(:,k-1)/self.V_g);
                    % X(:,k) = X(:,k-1) + Ts*(-p_2*X(:,k-1) + p_2*S_I*(I - I_b));
                    X(:,k) = self.A_dx*X(:,k-1) + self.B_dx*(I-self.I_b);
                    I_array(:,k) = self.A_di*I_array(:,k-1) + self.B_di*u_i(:,k-1);

                end
                x = [G; X; I_array];
                disp('1')
            end

        end
        function xhatm= dynamics(self,x,k)

            R_a = self.R_a(:,k-1);
            u_i = self.u_i(:,k-1);

            G = x(1);
            X_I = x(2);
            I_sc1 = x(3);
            I_sc2 = x(4);
            I_p = x(5);
            I_array = [I_sc1; I_sc2; I_p];

            G_hat = G + self.Ts*(-(self.S_g+X_I)*G + self.S_g*self.G_b + R_a/(self.V_g));
            % X_hat = X_I + Ts*(-p_2*X_I + p_2*S_I*(I_p/(BW*V_I) - I_b));
            % I_hat1 = I_sc1 + Ts*(-k_d*I_sc1 + u_i);
            % I_hat2 = I_sc2 + Ts*(-k_d*I_sc2 + k_d*I_sc1);
            % I_hatp = I_p + Ts*(-k_cl*I_p + k_d*I_sc2);
            X_hat = self.A_dx*X_I + self.B_dx*(I_p/(self.BW*self.V_I)-self.I_b);
            I_hata = self.A_di*I_array + self.B_di*u_i;

            % xhat = [G_hat; X_hat; I_hat1; I_hat2; I_hatp];
            xhatm = [G_hat; X_hat; I_hata];
        end
    end
end