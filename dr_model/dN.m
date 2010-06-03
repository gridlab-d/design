%% dN
% 
% Implements the aggregate demand response model for three types of DR
% control:
%
% delta - delta control affect the control band 0-L by moving it up (D>0) or
% down (D<0). The corresponds to thermostat setpoint reset.
%
% eta - direct load control affect the rate at which devices turn on
% (eta>0) or off (eta<0).
%
% phi - duty-cycle control affect the duty cycle of devices.
%
% The populations must be given as two arrays of the same length L with the
% number of devices in the off and on states at each temperature bin 0-L.
%
%%
function [dNoff dNon] = dN(Noff, Non, delta, eta, phi)

% check delta
if delta~=0
    error 'delta control not implemented yet';
end

% check eta
if eta < -1 || eta > 1
    error 'eta must be between -1 and 1'
end

% check phi
if phi <= 0
    error 'phi must be positive';
elseif phi <= 0.5 
    roff = phi/(1-phi);
    ron = 1;
elseif phi < 1
    roff = 1;
    ron = (1-phi)/phi;
else 
    error 'phi must be less than 1';
end

L = length(Non);
x = 1:(L-1);
if (eta>=0) % turning devices on
    dNon(L)    = -ron*Non(L)             + eta*Noff(L)   + (1-eta)*roff*Noff(L); 
    dNon(x)    = -ron*Non(x)             + eta*Noff(x)   + ron*Non(x+1);
    dNoff(1)   = -(1-eta)*roff*Noff(1)   - eta*Noff(1)   + ron*Non(1);
    dNoff(x+1) = -(1-eta)*roff*Noff(x+1) - eta*Noff(x+1) + (1-eta)*roff*Noff(x); 
else % (R<0) turning devices off
    dNon(L)    = -ron*(1+eta)*Non(L)     + eta*Non(L)    + roff*Noff(L);
    dNon(x)    = -ron*(1+eta)*Non(x)     + eta*Non(x)    + ron*(1+eta)*Non(x+1);
    dNoff(1)   = -roff*Noff(1)           - eta*Non(1)    + ron*(1+eta)*Non(1);
    dNoff(x+1) = -roff*Noff(x+1)         - eta*Non(x+1)  + roff*Noff(x);
end

