function [anp1,bnp1,historySum,tRBC_now] = get_historyContribution(t,an,bn,sn,zeroLoc,K,...
    dt,tRBC, W, snIdxNow, p,q,p0,q0,sn_hat,tRBC_now,phi_tn,phit_tn,s,tol)
% GET_HISTORYCONTRIBUTION a function to get the history contribution in BIE
% 
% INPUTs:
%  t: current time
%  an,bn: Fourier coefficients
%  sn,sn_hat: density grid function and transform
%  zeroLoc: location of the zero frequency
%  K: grid of frequency values 
%  dt: time-step
%  tRBC,tRBC_now: time-step for radiation boundary conditions and tRBC_now
%  is parameter to indicate which value on the time grid is being used at
%  the current step.
%  dt; time-step
%  W: width of the window function
%  snIdxNow: index of sn at the first time-step
%  p,q,p0,q0: the integrals p and q needed for the history evaluation 
%  (see notes). p0 and q0 correspond to the zero-frequency case. 
%  phi_tn,phit_tn: blending and window function for free-space projection
%  s:location of sources
%  tol: error tol
%  addProj: button to add free-space projection
%
% OUTPUTs: 
%  anp1,bnp1: Fourier coefficients at the new time-step
%  historySum: sum needed in the history evaluation on the right side
%  tRBC_now: the updated time to apply RBC

[anp1,bnp1] = alphaEvolution_spectral(an,bn,sn,zeroLoc,K,dt,W,snIdxNow,p,q,p0,q0,sn_hat);
if(t>=tRBC_now)
    [anp1,bnp1,~] = alphaRBCproj(anp1,bnp1,phi_tn,phit_tn);
    tRBC_now = tRBC_now + tRBC;
end

historySum = finufft1d2(s,-1,tol,anp1,struct('modeord',1));

end