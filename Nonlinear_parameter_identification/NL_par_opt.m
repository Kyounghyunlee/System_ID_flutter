clear all;
load ('ULCO_data_torsional_spring-1')

m=struct; % Creat a structure of measurement results;

p2=2.197186462*10^(-6); % Result from the normal form calculation (Unfolding parameter rescailing)
delta=U-flutter_speed;
nu=(-1+sqrt(1+4*p2*delta))/2/p2; % Coordinate transform for the unfolding parameter;

eig_vec=[0.0155 - 0.0160i;0.2474 + 0.2400i;-0.0000 - 0.0597i;0.9248 + 0.0000i;-0.0089 - 0.0034i;0.0528 - 0.1379i]; % Normalized vector of center space at Hopf bifurcation point
h_portion=abs(eig_vec(1))*2;
a=3; % Decide how many measurement set you are going to use. Note that measurement set=a+1 

%% Saving measurement results from the CBC
ml=length(U); 
m.nu=nu(ml-a:end); 
m.ff=flutter_freq(ml-a:end);
m.amp=amp(ml-a:end);
m.freq=freq(ml-a:end);
m.hportion=h_portion;
%%

fun = @(x)amp_LCO(x,m);
% fun = @(x)freq_LCO(x,m);
nonlcon = @constraint;
A = [];
b = [];
Aeq = [];
beq = [];
lb = [0,0]; % Lower bound of the parameter search
ub = [];
x0 = [3000,3000];
[x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon)

function [amp_err] = amp_LCO(x,m)
%Function that calculates the amplitude of the LCO
% x is the nonlinear stiffeness parameters
% nu is the unfolding parameter of hopf bifurcation
% amp_error = [amp_measured(ii)-amp_computed(ii)]^2/amp_measured(ii), where
% ii is the measurement set number

ml=length(m.nu);

% x(1): ka2, x(2)=ka(3)
for ii=1:ml
    r2=-0.00879561737287862*m.nu(ii)/(-7.10626449516451e-06*x(2) + 2.33467790333688e-08*x(1)^2);
    r=sqrt(r2);
    r=r*m.hportion;
    error(ii)=(r-m.amp(ii))^2/m.amp(ii)^2;
end

amp_err=sum(error);
end

function [c,ceq] = constraint(x)
%Nonlinear constraints for the optimization
%   x(1): ka2 , x(2): ka3
SLC=-7.10626449516451e-06*x(2) + 2.33467790333688e-08*x(1)^2; %Stability of the LCO near the hopf bifurcation point
% If SLC>0 LCO is unstable (subcritical Hopf bifurcation)
B23=8.57149841*10^(-10)*x(1)^2-2.615131149*10^(-6)*x(2); % If B23<0 frequency of LCO decreases as amplitude increase

c(1)= -SLC; % Find the parameters that makes dynamical system Subcritical Hopf Bifurcation (SLC>0)
ceq=[];
end


