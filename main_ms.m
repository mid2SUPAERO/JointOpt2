clear all; close all;
n_order=5;
PARAMETERS.F = 5000;                     % applied force in N
PARAMETERS.DT = 50;                      % variation of temperature
PARAMETERS.E1=210e3;                  % Young's modulus of adherend 1 in MPa
PARAMETERS.E2=70e3;                   % Young's modulus of adherend 2 in MPa
PARAMETERS.alpha1=12e-6;                % CTE of adherend 1 in K-1
PARAMETERS.alpha2=24e-6;                % CTE of adherend 1 in K-1
opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',2000,'MaxFunctionEvaluations',2e7);%'Display','iter',
% opts = optimoptions(@fmincon,'Algorithm','sqp');
% x0=[1000;0;-100;0;200];
x0=ones(n_order,1);

% problem = createOptimProblem('fmincon','objective',...
%  @(x)  my_func(x,n_order,PARAMETERS),'x0',x0,'lb',100*ones(1,n_order),'ub',2000*ones(1,n_order),'nonlcon',@(x)  non_linear_constraints(x,n_order,PARAMETERS),'options',opts);
% ms = MultiStart('UseParallel',true,'Display','iter');
problem = createOptimProblem('fmincon','objective',...
 @(x)  my_func(x,n_order,PARAMETERS),'x0',x0,'nonlcon',@(x)  non_linear_constraints(x,n_order,PARAMETERS),'options',opts);
ms = MultiStart('UseParallel',true,'Display','iter');
[x,f] = run(ms,problem,20);
[OBJECTIF]=post_processing(x,n_order,PARAMETERS);
switch PARAMETERS.E1
    case 210e3
        mat_string1='steel';
    case 70e3
        mat_string1='alu';
end
switch PARAMETERS.E2
    case 210e3
        mat_string2='steel';
    case 70e3
        mat_string2='alu';
end
Path=['order_',num2str(n_order),'_F_',num2str(PARAMETERS.F),'_DT_',num2str(PARAMETERS.DT),'_mat1_',mat_string1,'mat2',mat_string2,'/'];
mkdir(['order_',num2str(n_order),'_F_',num2str(PARAMETERS.F),'_DT_',num2str(PARAMETERS.DT),'_mat1_',mat_string1,'mat2',mat_string2])
save([Path,'results'],'x','f')
diary([Path,'diary.txt'])
print([Path,'results_plot'],'-dpng')


