%% dr_model_curtail
%
% Authors: DP Chassin, JC Fuller
% Written: May 19, 2010
% 
% This simulation illustrates a 10 MW curtailment in a 100 MW load
%

clear all; clc;

%% Simulation parameters
Qt      =   100;    % total load in MW
Qr      =    10;    % DR load in MW
eta     =  0.00;    % no initial demands/curtailments
DR      =  0.50;    % demand response signal (fraction of Qr)
phi     =  0.20;    % roff = 4 ron
dt      =     5;    % 5 minute/timestep;
db      =     2;    % 2F deadband control
P       =    45;    % cycle period in minutes
t0      =    24;    % time of curtailment (timestep)
t1      =    72;    % time of release (timestep)

%% Internal parameters
if phi < 0.5        % ron is higher rate 
    ron = 1;        % ron is unity rate
    roff = phi/(1+phi);
else                % roff is higher rate
    ron = (1+phi)/phi;
    roff = 1;       % roff is unity rate
end
L = ceil(P/dt/2);   % control band
T = 2*t1;   % timesteps to run

%% Initial state
Qu = Qt - Qr*phi;
Non  = ones(1,L)/L * phi;
Noff = ones(1,L)/L * (1-phi);

%% Output buffers
Q = zeros(T,1);
H = zeros(T,L);     % state histogram 
R = zeros(T,1);

%% Simulation run
for t=1:T
    
    % direct load control signal
    if t==1
        R(1) = eta;
    elseif t == t0
        R(t) = R(t-1)-DR;
    elseif t == t1
        R(t) = R(t-1)+DR;
    else
        R(t) = R(t-1);
    end
    
    % demand response
    [dNoff dNon] = dN(Noff, Non, 0, R(t), phi);
    
    % update population
    Noff = Noff + dNoff;
    Non = Non + dNon;
    
    % collect data
    H(t,:) = Non + Noff;
    Q(t) = Qu + Qr*sum(Non);
end

%% Output plots

figure(1);

subplot(2,1,1);
[AX,H1,H2] = plotyy((1:T)/(60/dt),Q,(1:T)/(60/dt),R);
xlabel('Time (h)');
set(get(AX(1),'Ylabel'),'String','Load (MW)');
set(get(AX(2),'Ylabel'),'String','Load control signal (\eta)');
title(sprintf('Period = %.0f min; Duty cycle = %.2f; Curtailment = %.0f%%', P, phi, DR*100));

% states
subplot(2,1,2);
surf(H); shading interp;
xlabel('State (x)');
ylabel('Time (t)');
zlabel('Population (n)');
title('Load states');
