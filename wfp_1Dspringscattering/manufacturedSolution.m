function [h0,dataParam] = manufacturedSolution(M,tFinal)
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

override = 0; % turn on '1' to avoid restrictions on chosen test function

mu = linspace(40,50,M)'; 
t0 = linspace(1,3,M)'; 

% choose t0 such that sig(0) = 0.01*eps
% t0 = sqrt(log(amp/(0.01*eps))./mu); 

%%% Working resolution using standard deviation 
N0 = 20;               % number of grid points per standard deviation (sd)
bumpWidth = sqrt(log(1/eps)./mu); 
h0_vec = bumpWidth./N0;
h0 = min(h0_vec); 

%%% Resolution computed from Fourier analysis
% omega0 = 2*sqrt(mu).*sqrt(log((1/tol)*sqrt(pi./mu))); 
% h0 = min(pi*(1 - gam)./omega0); 
% Nt0 = ceil(tFinal/h0); % or Nt0 = ceil(2*pi/h0);

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
