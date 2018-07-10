function [OBJECTIF]=my_func(x,n_order,PARAMETERS)
% global n_order
% =========================================================================
% 1D-BAR TEPS analysis
% single lap FGA joints 
% ==> in-plane tensile loading + uniform variation of temperature          
% ==> linear elastic materials
% author : Eric Paroissien
% date   : 23rd April 2018
% =========================================================================
% clear all
% close all
% clc
% =========================================================================
% =========================================================================
% =========================================================================
% INPUT DATA
% =========================================================================
% Size of the problem
% =========================================================================
n_K = 101;                    % total number of terms in the series (truncature at k=n_K-1)
% =========================================================================
% Geometry
% =========================================================================
ea = 0.2;                    % adhesive thickness in mm
e1 = 2.;                    % adherend 1 thickness in mm
e2 = 2.;                    % adherend 2 thickness in mm    
L = 25;                     % overlap length in mm
c = L/2;                    % half overlap in mm
b = 25;                     % overlap width in mm
% =========================================================================
% Material
% =========================================================================
% Adherends
E1 = PARAMETERS.E1;                  % Young's modulus of adherend 1 in MPa
E2 = PARAMETERS.E2;                   % Young's modulus of adherend 2 in MPa
alpha1 = PARAMETERS.alpha1;                % CTE of adherend 1 in K-1
alpha2 = PARAMETERS.alpha2;                % CTE of adherend 1 in K-1
% =========================================================================
% Adhesive
na_order=11;            % variable of each optimization case / has to be lower than n_K

Ga=zeros(n_K,1);        % variable to be optimized

% Ga(1,1)=x(1);
% Ga(3,1)=x(2);
% Ga(5,1)=x(3);
Ga(1:n_order)=x;

Kt=2;                   %

Ga_min=100;             % variable of each optimization case / has to be lower than n_K
Ga_max=2000;            % variable of each optimization case / has to be lower than n_K
contrainte_1=sum(Ga);   % has to be included between Ga_min and Ga_max
contrainte_2=Ga(1,1);   % has to be included between Ga_min and Ga_max
contrainte_3=Kt;
% =========================================================================
% Load
% =========================================================================
F = PARAMETERS.F;                     % applied force in N
DT = PARAMETERS.DT;                      % variation of temperature
% =========================================================================
% end of INPUT DATA
% =========================================================================
% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
% =========================================================================
% DEDUCED PARAMETERS
% =========================================================================
S1 = e1*b;                      % cross section area of adherend 1 (mm²)
S2 = e2*b;                      % cross section area of adherend 2 (mm²)
A1 = E1*S1;                     % membrane stiffness of adherend 1 (N)
A2 = E2*S2;                     % membrane stiffness of adherend 2 (N)
khi_A = A2/A1;                   % parameter for the adhrend membrane stiffness unbalance (-)
khi_alpha = alpha2/alpha1;       % parameter for the adhrend CTE unbalance  (-)
etaKRE = (1/e1/E1+1/e2/E2);     % governing parameter for the differential equation (mm/N)
kII=zeros(n_K,1);               % adhesive relative stiffness (N/mm3)
for i=1:na_order
    kII(i,1)=Ga(i,1)/ea;
end
% =========================================================================
% end of DEDUCED PARAMETERS
% =========================================================================
% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
% =========================================================================
% RECURSIVE MATRIX
% =========================================================================
A = zeros(n_K,n_K);
B = zeros(n_K,1);
for k=1:n_K-2
    for inc=1:k
        A(k,inc)=-c*c*etaKRE/(k-1+1)/(k-1+2)*kII(k-inc+1,1);
    end
    A(k,k+2)=1;
end
A(n_K-1,2)=1;
for inc=3:n_K
    A(n_K-1,inc)=inc-1;
end
A(n_K,2)=1;
for inc=3:n_K
    A(n_K,inc)=((-1)^(inc))*(inc-1);
end
B(n_K-1,1)=c*F/A2+c*(1-1/khi_alpha)*alpha2*DT;
B(n_K,1)=-c*khi_A*F/A2+c*(1-1/khi_alpha)*alpha2*DT;
% =========================================================================
% end of RECURSIVE MATRIX
% =========================================================================
% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
% =========================================================================
% COMPUTATION
% =========================================================================
U=A\B;
% =========================================================================
% end of COMPUTATION
% =========================================================================
% =========================================================================
% =========================================================================

% =========================================================================
% =========================================================================
% =========================================================================
% POST-PROCESSSING
% =========================================================================
% abscissa
% =========================================================================
m = 100;
X = zeros(2*m+1,1);
x = zeros(2*m+1,1);
X(1,1) = -1;
for i = 2:2*m+1
  X(i,1) = X(i-1,1)+1/(m);
end
for i = 1:2*m+1
  x(i,1) = X(i,1)*c;
end
Xik = zeros(2*m+1,n_K);
xik = zeros(2*m+1,n_K);
for i=1:2*m+1
    Xik(i,1)=1;
end
for i=1:2*m+1
    for k=2:n_K
        Xik(i,k)=X(i,1)^(k-1);
    end
end
for i=1:2*m+1
    xik(i,1)=1;
end
for i=1:2*m+1
    for k=2:n_K
        xik(i,k)=x(i,1)^(k-1);
    end
end
% =========================================================================
% results
% ========================================================================
u = zeros(n_K,1);
for i=1:n_K
    u(i,1)=U(i,1)/(c^(i-1));
end
G = zeros(2*m+1,1);
Du = zeros(m+1,1);                  % horizontal relative displacement in mm
T = zeros(m+1,1);                     % adhesive shear stress in MPa
G=Xik*Ga;
Du=xik*u;
for i = 1:2*m+1
    T(i,1) = G(i,1)/ea*Du(i,1);
end

OBJECTIF = max(T)-min(T);             % objectif à minimiser
Tmoy=5000/b/L;
contrainte_3=max(T)/Tmoy;

% =========================================================================
% Curves
% =========================================================================
% subplot(2,1,1)
% plot(x,G)
% legend('adhesive shear modulus')
% xlabel('abscissa along the overlap in mm')
% ylabel('modulus in MPa')
% grid on
% box on
% subplot(2,1,2)
% plot(x,T)
% legend('adhesive shear stress')
% xlabel('abscissa along the overlap in mm')
% ylabel('stress in MPa')
% grid on
% box on
% =========================================================================
% end of POST-PROCESSSING
% ========================================================================


