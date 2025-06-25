function phi = blending(phit, tau, tol)
% BLENDING retrieves the blending function
%
% phi = blending(phit, tau, tol) gives the value of blending function, the
%  antiderivative of the window function 'phit' at a given value 'tau' to
%  an error tolerance 'tol'

if nargin == 0, test_blending; return; end

phi = zeros(length(tau),1);
for i = 1:length(tau)
    % can replace quadgk by integral (not sure which is better)
    phi(i) = quadgk(phit,0, tau(i),'RelTol',tol,'AbsTol',tol); 
end

phi = real(phi);
end

function test_blending
tFinal = 1;
Nt = 100;
tol = 1e-12;
b = log(1/tol);
gam = 0.5;
w = ceil(2*b/(pi*gam));
t = linspace(0,tFinal,Nt); dt = t(2) - t(1);
phit = window(dt,b,w);
phi = blending(phit, t, tol);
plot(t,phi, 'LineWidth',2);
title('blending function \phi');
xlabel('x'); ylabel('y');
end