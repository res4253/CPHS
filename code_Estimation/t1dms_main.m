%ssogmm mode or inequality or linear
%modeはt_spanのシミュレーション時間になる(mode推定によるcutting_inds)
%inequalityの方はuva/padovaと同じシミュレーション時間でできる
% close all
clear;
clc;

Path = [
    './t1dms_data/subjects/:',...
    './t1dms_data/simu_data/:'];


Path = [
    '.\t1dms_data\subjects\;',...
    '.\t1dms_data\simu_data\;'];
addpath(Path);

%%%
d = "ssogmm"; %ssogmm inequality linear
for subject_num = 9
    for seed = 9


        if subject_num == 10
            load(['adult#0',num2str(subject_num),'.mat'],'BW','Ib','Gb');
        else
            load(['adult#00',num2str(subject_num),'.mat'],'BW','Ib','Gb');
        end

        load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed),'.mat'])

        % load("test008.mat")

        %%%
        Ts = 5;
        n_day = 4;
        t_end = 1440*n_day;
        N=size(modes,2);
        t_span = 200:Ts:200+N*5-1;

        %data procssing
        % modes = data_ssogmm.mode;
        if Ts==1
            modes = repelem(modes,5);
        elseif Ts==2.5
            modes = repelem(modes,2);
        end

        %%%
        G = data.G.signals.values';
        R = data.G.signals.values(1:Ts:end,:)';

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

        f=figure;
        plot(G,'LineWidth',2)
        hold on
        plot(data_ssogmm.ts,data_ssogmm.Gs,'LineWidth',2)
        legend('uva/padova','ssogmm')
        xlim([0 2600])
        grid on
        xticks(0:500:2500)   % 0, 2, 4, 6, 8, 10 に主目盛りを設定
yticks(100:50:250) % -1, -0.5, 0, 0.5, 1 に設定

        %%%%%
         
        if d == "inequality"
            %% mode %%
            Q = [data_ssogmm.xs(3,:);data_ssogmm.xs(4,:)];
            submodel = GastrointestinalSubmodelCL(p_esti,q_esti);

            Ra = [];
            inequality_modes = [];
            for t = 1:1153
                [Ra(t),inequality_modes(t)] = submodel.getR_a(Q(:,t),data_ssogmm.Ds(t));
            end

            f=figure;
            stairs(inequality_modes,'LineWidth',2)
            hold on
            stairs(cutting_inds,modes,'LineWidth',1)
            ylim([-0.5 1.5])
            legend('不等式','モード推定')
        end

        %%%                 %%%
        pass = './t1dms_data';
        filename = 'temp.mat';
        
        save(fullfile(pass,filename),'data_ssogmm',"subject_num","seed","Ts");

    end
end

rmpath(Path);