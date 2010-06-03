%% dr_price.m
%
% Generalized logit function to determine DR price
%
clear all;

eta = -0.1;
phi = 0.5;
phi = 1/(1-(phi-0.5)/0.5);
P = 50;

x = -100:0.01:100;

y = (1+phi*exp(-eta*(x-P))).^(-1/phi);

q = 0:0.01:1;

p = P-(log(q.^(-phi)-1) - log(phi))/eta;

plot(y,x,q,p,[0 1],[P P],':',[phi phi],[-100 100],':');
xlabel('Quantity');
ylabel('Price');
legend('Q = [1+\phi e^{-\eta(x-P_m)}]^{-1/\phi}','P = P_m - (log(Q^{-\phi}-1) - log(\phi))/\eta',3);
xlim([0 1]);
ylim([-100 200]);