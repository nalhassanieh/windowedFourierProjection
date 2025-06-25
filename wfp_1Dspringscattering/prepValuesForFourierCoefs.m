function [p,q,p0,q0,sn_hat,phi_tn,phit_tn] = prepValuesForFourierCoefs(sn,snIdxNow,dt,K,N,W,gl,s,tol,winData)
% PREPVALUESFORFOURIERCOEFS prepares parameters for the evaluation of the
% history part of the solution. 
%
% OUTPUT: 
%   the integrals p and q needed for the history evaluation (see notes). p0
%   and q0 correspond to the zero-frequency case.
%   sn_hat = \sum_{j = 1}^M \sigma_j(\tau)e^{ikx_j} at initial tau =
%   (0:W)*dt
%
% INPUT: 
%  'dt' timestep, 'K' vector of frequencies, 'W' width
%   of the window function, 'gl' GL nodes and weights on [-1,1], the
%   window function and its derivative 'phit' and 'phitt'; 's' is a vector
%   containing source locations.
%   'chebApproxInfo': a struct with fields domain [ax,bx] over which 
%      cheb weights are computed; ninters, number of subintervals, and nord, 
%      order of polynomial interpolation or number of cheb nodes   
%
% NOTE: 'expMat' and 'expMat2' 
%   correspond to matrices involving exponential evaluations needed in 
%   history treatment (needed in earlier version of the code)
  
% Get GL nodes and weights on the interval [0,dt]
[tau,w] = glwt(0,dt,gl);  

delta = W*dt; 

k = K';
tau = tau';
lK = length(K);

% allocate space for outputs
p = zeros(lK,(W+1)); q = zeros(lK,(W+1));
p0 = zeros(W+1,1); q0 = p0; 

for i = 0:W
    gam = i*dt;

    winArg = (tau + gam)./delta;
    [~,phit,phitt] = generalwindow(winArg,winData);
    phit = (1/delta)*phit;
    phitt = (1/delta^2)*phitt;

    I = (2*cos(k*(tau+gam)).*phit+(1./k).*(sin(k*(tau+gam)).*phitt));
    p(:,i+1) = (((1./k).*sin(k*(dt - tau))).*I)*w;
    q(:,i+1) = ((cos(k*(dt - tau))).*I)*w;

    I0 = (2*phit + (tau+gam).*phitt);
    p0(i+1) = sum(((dt - tau).*I0)'.*w);
    q0(i+1) = sum(I0'.*w);
end

sn_hat = finufft1d1(s,sn(snIdxNow-(0:W),:)',1,tol,N,struct('modeord',1));

%%% free-space parameters
phi_grid = (1:W) - 0.5;
[phi_tn1,phit_tn1,~] = generalwindow(phi_grid./W,winData);
phi_tn = phi_tn1'; phit_tn = (1/W)*phit_tn1';

end