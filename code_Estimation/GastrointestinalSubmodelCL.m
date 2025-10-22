classdef GastrointestinalSubmodelCL
    properties
        k_tau, k_abs0, k_abs1, BW, f, x_0, q_1, q_2
    end
    methods
        %これでいいのか？必要なのか？
        function self = GastrointestinalSubmodelCL(p,q_esti)
            self.k_tau = p(4);
            self.k_abs0 = p(5);
            self.k_abs1 = p(6);
            self.BW = p(12);
            self.f = p(13);
            self.x_0 = [0; 0];
            self.q_1 = q_esti(1);
            self.q_2 = q_esti(2);
        end

        function dx = dynamics(self,~,x,~,d)
            k_tau = self.k_tau;
            q_1 = self.q_1;
            q_2 = self.q_2;

            Q_1 = x(1);
            Q_2 = x(2);
            Q = Q_1 + Q_2;

            u_m = d(1);
            Ds = d(2);

            if Q > q_2*Ds || Q < q_1*Ds
                k_abs = self.k_abs1;
            else
                k_abs = self.k_abs0;
            end


            dQ_1 = -k_tau*Q_1 + u_m;
            dQ_2 = -k_abs*Q_2 + k_tau*Q_1;

            dx = [dQ_1; dQ_2];

        end

        function [R_a,mode] = getR_a(self, Q, Ds)
            BW = self.BW;
            f = self.f;
            q_1 = self.q_1;
            q_2 = self.q_2;

            Q_sum = Q(1)+Q(2);

            if Q_sum > q_2*Ds || Q_sum < q_1*Ds
                k_abs = self.k_abs1;
                mode = 1;
            else
                k_abs = self.k_abs0;
                mode = 0;
            end
            R_a = f .* k_abs .* Q(2) ./ BW;
        end

    end
end