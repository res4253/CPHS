close all
clear;
clc;

Path = [
    './t1dms_data/subjects/:',...
    './t1dms_data/simu_data/:'];

% Path = [
%     '.\t1dms_data\subjects\;',...
%     '.\t1dms_data\simu_data\;'];

addpath(Path);


%%%%

for subject_num =3 %subject
    for seed = 9 %seed

        d = "linear"; %ssogmm inequality linear


        if subject_num == 10
            load(['adult#0',num2str(subject_num),'.mat'],'BW','Ib','Gb');
        else
            load(['adult#00',num2str(subject_num),'.mat'],'BW','Ib','Gb');
        end

        load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed),'.mat'])

        %seed9によるパラメータ推定結果を利用
         % load(['t1dms_adult',num2str(subject_num),'_seed9.mat'],'p_esti','x_esti')


        %%%
        Ts = 5;
        n_day = 4;
        t_end = 1440*n_day;

        %data procssing
        % modes = data_ssogmm.mode;
        if Ts==1
            modes = repelem(modes,5);
        elseif Ts==2.5
            modes = repelem(modes,2);
        end

        N=size(modes,2);
        t_span = 200:Ts:200+1075*5-1;
        N_step = 1:N;

        %%%
        G = data.G.signals.values';

        insulin = data.injection.signals.values';
        meal = 1000*data.CHO.signals.values;

        %%%   ssogmm
        model = ssogmmCL("u_inb",(1/6)*Ib,...
            "BW",BW,...
            "G_b",Gb);

        model.p_0 = p_esti;
        model.x_0 = [x_esti(:,1)];%

        data_ssogmm = simu_t1dm_ssogmm(Ts,t_end,t_span,modes,...
            "insulin",insulin,...
            "meal",meal,...
            "ssogmm",model,"dynamics",d,"q",q_esti);

        %%%%%

        %discrete simulation
        %%%common setting
        Xs = data_ssogmm.xs;

        %estimate Xinit
        x_ini = x_esti(:,500);
        x_ini = x_ini([1 2 5 6 7]);

        %feedback
        Y = G(data_ssogmm.ts);
        % Y = data_ssogmm.Gs;

        %%input
        p = p_esti;
        u_m = data_ssogmm.meal;
        u_i = data_ssogmm.insulin;
        modes = data_ssogmm.modes;

        %%Ra
        R_a = get_Ra(Xs(4,:),p,modes,N);

        [xe,A_d,B_d1,B_d2,C,i_b,Ge] = linear_matrix(p,Ts);
        U = [u_m;u_i;R_a;modes];
        

        %%%

        %%%linear
        delta_x=zeros(5,N);
        x_linear=zeros(5,N);

        delta_x(:,1) = Xs([1 2 5 6 7],1)-xe;
        x_linear(:,1) = Xs([1 2 5 6 7],1);

        for k=2:N

            delta_u_i = u_i(:,k-1)-i_b;

            delta_x(:,k) = A_d*delta_x(:,k-1) + B_d1*delta_u_i + B_d2*R_a(:,k-1);
            x_linear(:,k) = delta_x(:,k) + xe;
        end


        %observer
        % 状態推定誤差の重み
        Q = diag([1, 1, 1, 1, 1]);%
        R = 1e16;
        [L,pole] = gain(p,Ts,Q,R);

        delta_xob = zeros(5,N);
        xob = zeros(5,N);

        delta_xob(:,1) = x_ini - xe;
        xob(:,1) = x_ini;

        for k=2:N

            delta_y = Y(:,k-1) - Ge;
            delta_u_i = u_i(:,k-1) - i_b;

            delta_xob(:,k) = A_d*delta_xob(:,k-1) + B_d1*delta_u_i + B_d2*R_a(:,k-1) + L*(delta_y - C*delta_xob(:,k-1));
            xob(:,k) = delta_xob(:,k) + xe;
        end



        %%ssogmm
        dic = dic_ssogmmCL(p,Ts,N,U);
        x_ssogmm = dic.dicrete(Xs,"ZOH");


        %kalman
        A = dic.A;
        A_di = dic.A_di;
        A_dx = dic.A_dx;
        C = [1 0 0 0 0];

        xhat = zeros(5,N);

        a = [1,1,1,1,1];
        P = 1e-6*diag(a);
        xhat(:,1) = x_ini;

        %%
        Q = eps; %システム雑音
        R = eps; %観測雑音
        % v = randn(1,N)*sqrtm(Q);
        % w = randn(1,N)*sqrtm(R);
        B = [1,1,1,1,1];
        B_diag = diag(B);
        g_block = [];
        v = [];
        v = [v,Y(1,1)-x_ini(1,1)];
        for k = 2:N

            xhatm = dic.dynamics(xhat(:,k-1),k);
            v = [v,Y(:,k) - xhatm(1,1)];

            A_tmp = [zeros(3,2) A_di];
            A_temp =[A(xhat(1,k-1),xhat(2,k-1)) ; A_tmp];
            % Pm = A(xhat(1,k-1),xhat(2,k-1))*P*A(xhat(1,k-1),xhat(2,k-1)) + b*Q*b';
            Pm = A_temp*P*A_temp' + B_diag*Q*B_diag';

            g = Pm*C'/(C*Pm*C' + R);

            xhat(:,k) = xhatm + g*(Y(:,k) - C*xhatm);
            P = (eye(5) - g*C)*Pm;
            g_block = [g_block,g];
        end
 
    end
end


ssogmm = data_ssogmm.xs([1 2 5 6 7],:);
for i=1:5
    figure
    plot(ssogmm(i,:),'LineWidth',2)
    grid on
    ax = gca;
    saveFig(ax,'png')
end