%% dr_model_pid_phi
%
% cyclic model of a thermostat cooling device supporting regulation with
% PID control of the rate phi
%

clear all; clc;

% load model
T =    100;         % number of timesteps
N =    100;         % number of loads (thousands)
q =      5;         % load capacity (kW)
phi = 0.20;         % load natural duty cycle )0..1(
eta = 0.01;         % demand probability (/t)
L =     10;         % load control band (0=off, L=on)

% control model
Pu = 10;             % critical damping period (zero disables controls)
if (Pu>0)
    Kp = 0.12;          % proportional gain
    Ki = 2.0*Kp/Pu;     % integral gain
    Kd = Kp*Pu/8;       % differential gain (zero to disable)
else 
    Kp = 0;
    Ki = 0;
    Kd = 0;
end

% other system fluctuations (net)
V =   0.2;          % random fluctuations (MW)
M =   0.5;          % random walk (MW/t)

% disturbance
Ts =   10;          % start time of disturbance (t, 0=none) 
Se =  -10;          % size of disturbance (MW)
Te =   30;          % end time of disturbance (t, 0=none)


%% initial conditions
if phi < 0
    error "phi must be positive";
elseif phi <= 0.5 
    roff = phi/(1-phi);
    ron = 1;
elseif phi <= 1
    roff = 1;
    ron = (1-phi)/phi;
else 
    error "phi must be less than or equal to 1";
end

% initial control state
R = phi;            % initial control is 0
S = phi*N*q;        % scheduled load (remains constant)

% initial state populations
x = 0.5:(L-0.5);      % bins
if (eta==0) % uniform distribution
    Non = N*ones(1,L)*phi/L;
else % exponential distribution
    Non = N * eta * (1-phi) * exp(eta*(L-x)/roff) / (exp(eta*L/roff)-1);
end
Noff = Non * ron / roff;

%% output buffers
G = zeros(T,1);     % generation
Q = zeros(T,1);     % total load
D = zeros(T,1);     % dispatch response 
H = zeros(T,L);     % state histogram 
W = zeros(T,1);     % population status

% controller inputs
E = zeros(T,1);     % imbalance
Ei = 0;             % imbalance integral
Ed = 0;             % imbalance differential

%% simulation
for t=1:T
    
    % events
    if (t==Ts)
        S = S + Se;
    elseif (t==Te)
        S = S - Se;
    end
    S = S + randn*M;
    G(t) = S + randn*V;
    
    % total load calc
    Q(t) = sum(Non)*q;    % total load
    
    % imbalance calc
    E(t) = S - Q(t);      % imbalance 
    
    % eta is controllable based on imbalance as a fraction of schedule
    Ep = E(t)/S;
    Ei = Ei + Ep;
    if t>1
        Ed = (E(t) - E(t-1))/S;
    end
    
    % PID control
    R = phi + Kp*Ep + Ki*Ei + Kd*Ed; 
    
    % limit saturation
    if R < 0 
        R = 0; 
        disp(sprintf('duty-cycle control saturation (t=%.0f, R=0)',t));
    elseif R > 1 
        R = 1; 
        disp(sprintf('duty-cycle control saturation (t=%.0f, R=+1)',t));
    end

    % collect waterfall histogram data or dispatch history
    H(t,:) = Non + Noff;
    D(t) = (R - phi);
    
    % device state transitions
    [dNoff dNon] = dN(Noff, Non, 0, eta, R);

    % population update
    Non = Non + dNon;
    Noff = Noff + dNoff;
    
    % population status
    W(t) = sum(Non)/N;
    
end

%% output
figure(2);

% generation profile
subplot(2,2,1);
plot(1:T,G,1:T,Q);
xlabel('Time (t)');
ylabel('Power (MW)');
title('Power');
legend('Generation','Load',0);

% states
subplot(2,2,2);
surf(H); shading interp;
xlabel('State (x)');
ylabel('Time (t)');
zlabel('Population (n)');
title('Load states');

subplot(2,2,3);
[AX,H1,H2] = plotyy(1:T,E./G*100,1:T,D);
xlabel('Time (t)');
set(get(AX(1),'Ylabel'),'String','Imbalance (%)');
set(get(AX(2),'Ylabel'),'String','Duty-cycle control signal (\phi)');
title('Feedback control');

% population status
subplot(2,2,4);
plot(W);
xlabel('Time (t)');
ylabel('Fraction of active load (%)');
title('Load status');

disp('Done.');
