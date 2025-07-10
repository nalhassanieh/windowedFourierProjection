function [dt0,dataParam] = get_uIncidentInfo(dataParam)
% GET_UINCIDENTINFO prepares information needed for an incident wave of the
% form uin = exp(-mu*(x - t - t0).^2);
%
% INPUT: 
%  tFinal: final time
%  tol: error tolerance for the full problem (to compute starting dt)
%  gam: pad gam*Nyquist for window function (to compute starting dt)
%
% OUTPUT: 
%   mu,t0 as indicated above
%   Nt0: number of initial time steps to resolve incident pulse 
%   h0:  minimum step size to resolve incident pulse
%   dataParam: is struct to hold mu and t0 (for ease of switching between
%   types of test solutions

if nargin == 0, test_get_uIncidentInfo(); return; end

mu = dataParam.mu; 
t0 = dataParam.t0;

if(isnan(t0))
    t0 = -sqrt(log(1/eps)./mu);
end

% choose the initial time step
dt0 = pi/(4*sqrt(mu*log(1/eps)));

if(isfield(dataParam,'doubleTimeStep'))
    if(dataParam.doubleTimeStep == 1)
        dt0 = 2*dt0;
    end
end

dataParam.t0 = t0;

end

function test_get_uIncidentInfo()

mu = 5; 
[dt0,dataParam] = get_uIncidentInfo(mu); 

x = linspace(-pi,pi,1000); t = linspace(0,3*pi,1000); t0 = dataParam.t0; 
uin = get_uIncident(x,t,mu,t0); 

figure(1)
imagesc(x,t,uin);

end