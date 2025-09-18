% % MainLoop function for WSAO and draw curve 
% % one run for all CEC2017 data and 30 repeat_times 
% % It takes much time like 10-30 minutes for all CEC2017 benchmark
% function,but You can stop the program at any time without causing serious consequences.
clear;
clc;

TEST_fun = [1,3:30]';
all_avebest = zeros(30,1);
all_std = zeros(30,1);
for fun_num = [1,3:30]
    Particles_no=100 ; % Number of search agents
    
    %Function_name='F12';
    Function_name = ['F',num2str(fun_num)];         % CEC2017(F1、F3~F30) 
    Max_iteration=500; % Maximum numbef of iterations
    Min_DomainR=0;
    
    repeat_times = 30;
    totol_Best_score=[];
    totol_Best_pos=[];
    totol_Curve=[];
    
    % Load details of the selected benchmark function
    [lb,ub,dim,fobj]=CEC2017(Function_name);
    
    for i = 1:repeat_times
    
        [Best_score,Best_pos,WestPSO_cg_curve]=WSAO(Particles_no,Max_iteration,lb,ub,dim,fobj);
        totol_Best_score=[totol_Best_score;Best_score];
        totol_Best_pos=[totol_Best_pos;Best_pos];
        totol_Curve=[totol_Curve;WestPSO_cg_curve];
    end

    %保存数据
    mkdir("data_save")
    save(['./data_save/F',num2str(fun_num)],'totol_Curve')
    all_avebest(fun_num) = mean(totol_Best_score);
    all_std(fun_num) = std(totol_Best_score);
end

save('./data_save/all_avebest','all_avebest');
save('./data_save/all_std','all_std');
