function [dt0,dataParam] = manufacturedSolution(M,tFinal)
% MANUFACTUREDSOLUTION  analytic density functions
%
% [sig,Nx0,h0,t0,mu] = manufacturedSolution(M,tFinal,s,beta)
%  returns in sig a cell array of M density function handles
%  
% Input values: 
%  M:      number of spring scatterers (sources)
%  tFinal: final time
%  tol: error tolerance for the full problem (to compute starting dt)
%  gam: pad gam*Nyquist for window function (to compute starting dt)
%
% Outputs:
%  sig: function handle in terms of time that produces a vector of density
%       values evaluated at a specific time t: 
%       [sig_1(t),...,sig_M(t)]
%  Nx0, h0: Number of grid points needed to ensure sig is resolved on the
%           time grid, as well as the maximum step-size.
%  t0, mu : (to document)
%  dataParam: is a struct involving t0,mu and sig.
%
%  Notes:
% Take the densities to be Gaussians for testing. 

if nargin == 0, test_manufacturedSolution(); return; end

override = 0; % turn on '1' to avoid restrictions on chosen test function

% mu = linspace(40,50,M)'; 
mu = linspace(2,5,M)'; 

% choose t0 such that sig(0) = eps
t0 = sqrt(log(1/eps)./mu); 

% choose the initial time step
dt0_vec = pi./(4*sqrt(mu*log(1/eps)));
dt0 = min(dt0_vec); 

% ensure that the forcing function is bandlimited within the computational time domain
if(any(2*t0>tFinal)&&(override == 0))
    error('Forcing inappropriate for testing. For tFinal=%1.2f\nChoose other mu or tFinal',tFinal);
end

% Define the density functions to be Gaussians
sig = @(t) exp(-mu.*((t - t0).^2));
   
dataParam.sig = sig; 
dataParam.t0 = t0; 
dataParam.mu = mu; 

end

%%%% 
function test_manufacturedSolution()

M = 1; tFinal = 3*pi; 
[dt0,dataParam] = manufacturedSolution(M,tFinal); 

t = linspace(0,tFinal,1000);
figure(1)
plot(t,dataParam.sig(t)); grid on; 

end