%やらなくていいかも
%クビです


% close all
clear;
clc;

Path = ['/home/tanabe/code/miwa/model/:',...
    '/home/tanabe/code/miwa/database/:',...
    '/home/tanabe/code/t1dms_data/:',...
    '/home/tanabe/code/t1dms_data/subjects/:',...
    '/home/tanabe/code/t1dms_data/simu_data/:'];
addpath(Path);

%%%

for j = 1
seed = 1;
subject_num = j;

if subject_num == 10
    load(['adult#0',num2str(subject_num),'.mat'],'BW');
    load(['t1dms_adult0',num2str(subject_num),'_seed',num2str(seed),'.mat'])
else
    load(['adult#00',num2str(subject_num),'.mat'],'BW');
    load(['t1dms_adult00',num2str(subject_num),'_seed',num2str(seed),'.mat'])
end


%%%
Ts = 5;
n_day = 4;
t_end = 1440*n_day;
N=1075;
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
R = data.G.signals.values(1:Ts:end)';
cutting_inds = (cutting_inds(1)+1):(cutting_inds(end)+1);
R = R(cutting_inds);

insulin = data.injection.signals.values';
meal = 1000*data.CHO.signals.values;

%%%   ssogmm
model = ssogmmCL("u_inb",(1/6)*insulin(1,1),...
    "BW",BW,...
    "G_b",G(1,1));

model.p_0 = p_esti;
model.x_0 = [100;0;x_esti(3,1);x_esti(4,1);5000;1500;300];%

L = gain(model.p_0,0);

data_obs = simu_obs(Ts,t_end,t_span,modes,R,L,...
                    "insulin",insulin,"meal",meal,"ssogmm",model);

figure
plot(G,'LineWidth',2)
hold on
plot(data_obs.ts,data_obs.Gs,'LineWidth',2)


%rmse all
xob = data_obs.xs;
Error_all = rmse(xob(1,:)',R');

%rmse only 4day
N_step =1:N;
d = floor(size(N_step,2)/4)*3;
Error_4day = rmse(xob(1,d:end)',R(:,d:end)');

rmpath(Path);


end