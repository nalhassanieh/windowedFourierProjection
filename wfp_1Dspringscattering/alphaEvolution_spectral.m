function [anp1,bnp1] = alphaEvolution_spectral(an,bn,sn,zeroLoc,K,dt,W,snIdxNow,p,q,p0,q0,sn_hat)
% ALPHAEVOLUTION_SPECTRAL computes the Fourier coefficients to spectral
% accuracy
%
% [anp1,bnp1] = alphaEvolution_spectral(an,bn,sn,zeroLoc,K,dt,W,snIdxNow,...
%               p,q,p0,q0,s,N,tol)
% returns Fourier coefficients at the new time t = n*dt, anp1,...
% and time derivativesof anp1, bnp1. 
% 
% Inputs: 
%  snIdxNow : the current time index
%  an: Fourier coefficients at time t = tInitial + (n-1)*dt
%  bn: 1st time derivative of an
%  sn: matrix of density grid functions sn(time, source number)
%  zeroLoc: location of the zero wave number k
%  K: vector of fourier wave numbers
%  N: number of Fourier modes
%  tol: error tolerance for finufft
%  dt: time-step
%  W: W*dt = delta where [0,delta] is the support of phit 
%  s: [s(1),...,s(M)] location of the M spring scatterers
%  p,q,p0,q0: integrals needed in the evaluation of Fourier coefficients
%  (see notes)
%  expMat: expMat = exp(1i*k*(s')) needed in Fourier coefficient evaluation

if nargin == 0, test_alphaEvaluation_spectral; return; end

%% Treatment of u_history
k = K';

% update k values not equal to zero
i = 0:W; 
h = sum(p(:,i+1).*sn_hat,2); % sum rows
g = sum(q(:,i+1).*sn_hat,2);
h = h*dt; g = g*dt;

% evolution equations k ~= 0
anp1 = an.*cos(k*dt) + bn.*(1./k).*sin(k*dt) + h;
bnp1 = -k.*an.*sin(k*dt) + bn.*cos(k*dt) + g;

% update k = 0 frequency 
sig0 = sum(sn(snIdxNow-i,:),2); % sum rows
h0 = sum(p0(i+1).*sig0);
g0 = sum(q0(i+1).*sig0);
h0 = h0*dt; g0 = g0*dt;

% evolution equations for k = 0 case
anp1(zeroLoc) = an(zeroLoc) + bn(zeroLoc)*dt + h0;
bnp1(zeroLoc) = bn(zeroLoc) + g0;

end