function [an,bn,sn,sn_hat] = prepForNextStep(anp1,bnp1,sn,sn_hat,snIdxNow,tgh,s,N,tol)
% PREPFORNEXTSTEP prepares for the next step in the the time-stepping
% 
% INPUTs: 
%  anp1,bnp1: Fourier coefficients at new step and time derivative
%  sn, sn_hat: density grid function and transfrom
%  snIdxNow: index of density at the first time step
%  tgh: number of ghost points in time
%  s: location of sources
%  N: number of Fourier coefficients
%  tol: error tolerance
%
% OUTPUTs:
%  an,bn: Fourier coefficients at the old step
%  sn,sn_hat: density grid function and transfrom

% update values for the next time step
an = anp1; bn = bnp1;

% Update the density values for the next time-step
for kk = 1:(tgh + 1)
    sn(kk,:) = sn((kk+1),:);
end
sn_hat1 = finufft1d1(s,sn(snIdxNow,:)',1,tol,N,struct('modeord',1));
sn_hat = [sn_hat1,sn_hat(:,1:(end-1))];

end