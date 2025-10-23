close all
clear
clc


for i = 51:60
    adult(i).q = [];
    adult(i).p = [];
    for j = 1:10

        subject_num = i;
        seed = j;

        % load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed)])
        load(['id',num2str(subject_num),'_seed',num2str(seed)])


        adult(i).q = [adult(i).q,q_esti];

        adult(i).p = [adult(i).p,p_esti];

        
    end
    adult(i).mean = mean(adult(i).p,2);
    adult(i).median = median(adult(i).p,2);
    adult(i).std = std(adult(i).p,[],2);
    adult(i).Normalization = normalize(adult(i).p,2);

    % 
    % f = figure;
    % f.Position = [1271,450,560,420];
    % boxchart(adult(i).Normalization(:,:)')
    % title(['num,',num2str(subject_num)])

    f = figure;
    f.Position = [12,441,622,420];
    tiledlayout(2,5)
    for k = 2:11
        nexttile
        plot(adult(i).p(k,:))
        title(k)
        xline(6)
    end
end

% 


%ssogmm paramator
% % G_b = p(1);
% % V_I = p(2);%%
% % S_I = p(3);%
% % k_tau = p(4);
% % k_abs = p(5);
% % k_abs = p(6);
% % k_d = p(7);
% % k_cl = p(8);%%
% % S_g = p(9);
% % V_g = p(10);
% % p_2 = p(11);
% % BW = p(12);
% % f_c = p(13);
% % I_b = p(14);





% for i = 1:10
%     for j = 1:10
% 
%         subject_num = i;
%         seed = j;
% 
%         load(['t1dms_adult',num2str(subject_num),'_seed',num2str(seed)])
% 
%         f = figure;
%         f.Position = [28,101,1832,872];
%         tiledlayout(2,4)
%         for k = 1:7
%             nexttile
%             plot(x_esti(k,:))
%             title(k)
%         end
%         xlabel(['num',num2str(subject_num),'seed',num2str(seed)])
%     end
% end