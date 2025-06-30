function [dt0,dataParam] = get_uIncidentInfo(mu)
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

t0 = -3; 

% choose the initial time step
dt0 = pi/(4*sqrt(mu*log(1/eps)));

dataParam.mu = mu; 
dataParam.t0 = t0; 

end