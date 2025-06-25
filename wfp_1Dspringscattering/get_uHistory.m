function uh = get_uHistory(an,x,tol)
% GET_UHISTORY computes the history part of the solution 
%
% uh = get_uHistory(an,x,K) 
%  computes the history part of the numerical solution, given Fourier
%  modes 'an', 'x' grid points

%%% Using Finufft (NU to NU)
% uh = (1/(2*pi))*finufft1d3(K,an,-1,tol,x);

%%% Using Finufft (U to NU)
opts.modeord = 1; 
uh = (1/(2*pi))*finufft1d2(x,-1,tol,an,opts);

% Did not use fft because number of test points is small. We are using here
% N frequencies to compute Nx values where Nx<N. The fft operator truncates
% at Nx frequencies, which produces incorrect results. 

%%% Standard if missing Finufft
% uh = zeros(size(x));
% cnt = 1;
% for k = K
%     uh = uh + an(cnt)*exp(-1i*k*x);
%     cnt = cnt + 1;
% end
% uh = (1/(2*pi))*uh;

end