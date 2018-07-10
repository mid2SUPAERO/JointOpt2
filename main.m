%% Main
n_order=5;
PARAMETERS.F = 5000;                     % applied force in N
PARAMETERS.DT = 0;                      % variation of temperature
PARAMETERS.E1=210e3;                  % Young's modulus of adherend 1 in MPa
PARAMETERS.E2=210e3;                   % Young's modulus of adherend 2 in MPa
PARAMETERS.alpha1=12e-6;                % CTE of adherend 1 in K-1
PARAMETERS.alpha2=12e-6;                % CTE of adherend 1 in K-1
opts = optimoptions('fmincon','Algorithm','sqp','MaxIterations',2000,'MaxFunctionEvaluations',2e7);%'Display','iter',
% opts = optimoptions(@fmincon,'Algorithm','sqp');
% x0=[1000;0;-100;0;200];
x0=ones(n_order,1);
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxIterations',2000,'MaxFunctionEvaluations',2e7);
A = [];
b = [];
Aeq = [];
beq = [];

% lb = [ 100; 300 ];
% ub = [2000;6000];
lb = [];
ub = [];
x = fmincon(@(x) my_func(x,n_order,PARAMETERS),x0,A,b,Aeq,beq,lb,ub,@(x) non_linear_constraints(x,n_order,PARAMETERS),options);
[OBJECTIF]=post_processing(x,n_order,PARAMETERS);