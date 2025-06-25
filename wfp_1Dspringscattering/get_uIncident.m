function uin = get_uIncident(x,t,mu,t0)
% GET_UINCIDENT returns an incident wave
%
% Inputs: 
%   x,t: space and time variables
%   mu,t0: for the Gaussian as seen in formula below

if nargin == 0, test_get_uIncident; return; end  

uin = exp(-mu.*(x - t - t0).^2); 

end

function test_get_uIncident
ax = -pi; bx = pi; Nx = 100; 
x = linspace(ax,bx,Nx); 
t = 3; 
mu = 30; 
t0 = -3; 
uin = get_uIncident(x,t,mu,t0); 
plot(x,uin,'LineWidth',2); 
xlabel x; ylabel u; 
title 'plot of u_in'; 
end
