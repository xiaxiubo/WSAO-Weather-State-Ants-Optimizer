% % Main function of WSAO and draw curve 
clear;
clc;

Particles_no=100 ; % Number of search agents

Function_name='F15'; % Name of the test function that can be F1，F3 to F30 for CEC2017 (Table 2 in the paper)

Max_iteration=500; % Maximum numbef of iterations


% Load details of the selected benchmark function
[lb,ub,dim,fobj]=CEC2017(Function_name);

[Best_score,Best_pos,WestPSO_cg_curve]=WSAO(Particles_no,Max_iteration,lb,ub,dim,fobj);


figure('Position',[269   240   800   290])
%Draw search space
subplot(1,2,1);
func_plot_cec17(Function_name);
title('Parameter space')
xlabel('x_1');
ylabel('x_2');
zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
subplot(1,2,2);
semilogy(WestPSO_cg_curve,'Color','r')
% semilogy(totol_Curve','Color','r')
title('WSAO Objective space')
xlabel('Iteration');
ylabel('Best score obtained so far');

grid on
box on
legend('WSAO')

display(['The best solution obtained by WestPSO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by WestPSO is : ', num2str(Best_score)]);
