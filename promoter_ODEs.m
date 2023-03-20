function dydt = promoter_ODEs(t,y,t_trace,Msn2_trace,promoter_params,K_scale,fraction_active)

% Global Parameters
d3 = 0.08; % Turn off if modeling mRNA decay
k4 = 15;
% d4 = 0.003;
d4 = 0.001;
k5 = 0.06;

% Promoter specific parameters
k1 = promoter_params(1);
d1 = promoter_params(2);
k2 = promoter_params(3);
K = promoter_params(4);
n = promoter_params(5);
d2 = promoter_params(6);
k3 = promoter_params(7);

% ODE indices
dydt = zeros(6,1);
P_unbound = y(1);
P_bound = y(2);
P_active = y(3);
mRNA = y(4);
YFP = y(5);
mYFP = y(6);

% Get Msn2 value
Msn2 = interp1(t_trace,Msn2_trace,t);
% Msn2 = qinterp1(t_trace,Msn2_trace,t);

% Scale by fraction active
Msn2 = fraction_active*Msn2;

% 1 Hill and 1 Msn2 (scale K)
K = K_scale*K;
dydt(1) = d1*P_bound - ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound;
dydt(2) = ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound + d2*P_active - d1*P_bound - k2*Msn2*P_bound;
dydt(3) = k2*Msn2*P_bound - d2*P_active;
dydt(4) = k3*P_active - d3*mRNA;
dydt(5) = k4*mRNA - d4*YFP - k5*YFP;
dydt(6) = k5*YFP - d4*mYFP;

% 1 Hill and 1 Msn2 scale Kd (where Kd = K^n)
% Kd = K_scale*(K^n);
% dydt(1) = d1*P_bound - ((k1*Msn2^n)/(Kd + Msn2^n))*P_unbound;
% dydt(2) = ((k1*Msn2^n)/(Kd + Msn2^n))*P_unbound + d2*P_active - d1*P_bound - k2*Msn2*P_bound;
% dydt(3) = k2*Msn2*P_bound - d2*P_active;
% dydt(4) = k3*P_active - d3*mRNA;
% dydt(5) = k4*mRNA - d4*YFP - k5*YFP;
% dydt(6) = k5*YFP - d4*mYFP;

% 1 Msn2 and 1 Hill 
% dydt(1) = d1*P_bound - k1*Msn2*P_unbound;
% dydt(2) = k1*Msn2*P_unbound + d2*P_active - d1*P_bound - ((k2*Msn2^n)/(K^n + Msn2^n))*P_bound;
% dydt(3) = ((k2*Msn2^n)/(K^n + Msn2^n))*P_bound - d2*P_active;
% dydt(4) = k3*P_active - d3*mRNA;
% dydt(5) = k4*mRNA - d4*YFP - k5*YFP;
% dydt(6) = k5*YFP - d4*mYFP;

% 1 Hill and 2 Msn2
% dydt(1) = d1*P_bound - ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound;
% dydt(2) = ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound + d2*P_active - d1*P_bound - k2*Msn2*P_bound;
% dydt(3) = k2*Msn2*P_bound - d2*P_active;
% dydt(4) = k3*Msn2*P_active - d3*mRNA;
% dydt(5) = k4*mRNA - d4*YFP - k5*YFP;
% dydt(6) = k5*YFP - d4*mYFP;

% 2 Hill
% dydt(1) = d1*P_bound - ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound;
% dydt(2) = ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound + d2*P_active - d1*P_bound - ((k2*Msn2^n)/(K^n + Msn2^n))*P_bound;
% dydt(3) = ((k2*Msn2^n)/(K^n + Msn2^n))*P_bound - d2*P_active;
% dydt(4) = k3*P_active - d3*mRNA;
% dydt(5) = k4*mRNA - d4*YFP - k5*YFP;
% dydt(6) = k5*YFP - d4*mYFP;

% 3 Hill
% dydt(1) = d1*P_bound - ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound;
% dydt(2) = ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound + d2*P_active - d1*P_bound - ((k2*Msn2^n)/(K^n + Msn2^n))*P_bound;
% dydt(3) = ((k2*Msn2^n)/(K^n + Msn2^n))*P_bound - d2*P_active;
% dydt(4) = ((k3*Msn2^n)/(K^n + Msn2^n))*P_active - d3*mRNA;
% dydt(5) = k4*mRNA - d4*YFP - k5*YFP;
% dydt(6) = k5*YFP - d4*mYFP;

% Lee 2021 ODEs
% dydt(1) = d1*P_bound - ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound;
% dydt(2) = ((k1*Msn2^n)/(K^n + Msn2^n))*P_unbound + d2*P_active - d1*P_bound - k2*P_bound;
% dydt(3) = k2*P_bound - d2*P_active;
% dydt(4) = k3*P_active - d3*mRNA;
% dydt(5) = k4*mRNA - d4*YFP - k5*YFP;
% dydt(6) = k5*YFP - d4*mYFP;

% Hansen 2013 
% Msn2_max = 1;
% dydt(1) = d1.*P_bound - k1*(Msn2/Msn2_max)*P_unbound;
% dydt(2) = k1*(Msn2/Msn2_max)*P_unbound + d2*P_active - (d1 + (k2*Msn2^n)/(K^n + Msn2^n))*P_bound;
% dydt(3) = (k2*Msn2^n/(K^n + Msn2^n))*P_bound - d2*P_active;
% dydt(4) = (k3*Msn2^n/(K^n + Msn2^n))*P_active - d3*mRNA;
% dydt(5) = k4*mRNA - (d4 +k5)*YFP;
% dydt(6) = k5*YFP - d4*mYFP;

