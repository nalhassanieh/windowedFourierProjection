function y = chebEval(x,wts,chebApproxInfo)
% CHEBEval retruns the approximation of a function at x
%
% y = chebApprox(x,ax,h,wts,ninters,nord) approximates a function f
% (included in the wts argument) using chebyshev polynomial expansion on
% subintervals of a domain that starts at ax
%
% INPUTS:
%    x: value in the domain where f is to be approximated
%    wts: chebyshev weights multiplied by function values in all
%         subintervals, wts(nord,ninters)
%    chebApproxInfo.domain: domain [ax,bx] overwhich wts are computed. 
%    chebApproxInfo.ninters: number of subintervals
%    chebApproxInfo.nord: order of polynomial interpolation; number of cheb 
%         nodes
% NOTE: code adjusted from code provided by Leslie Greengard. 

if nargin == 0, test_chebEval; return; end

domain = chebApproxInfo.domain; ax = domain(1); bx = domain(2); 
ninters = chebApproxInfo.ninters; 
nord = chebApproxInfo.nord; 

Nx = length(x);
y = zeros(Nx,1); 

h = (bx-ax)/ninters; 

for ix = 1:Nx
    % get index of left node in current interval
    iint = min((floor((x(ix)-ax)/h) + 1),ninters);

    % get interval endpoints [a,b]
    a = ax + (iint-1)*h; b = ax + iint*h;

    % transform to u in [-1,1]
    u = (2*x(ix) - a - b)/(b-a);

    % Evaluate the cheb polynomial
    y(ix) = wts(1,iint);
    for J = 2:nord
        TJ = cos((J - 1) * acos(u));  % Chebyshev polynomial T_J(X)
        y(ix) = y(ix) + TJ * wts(J,iint); % Accumulate the value
    end
end
end

function test_chebEval
ninters = 12; % number of subintervals
nord = 4;     % number of cheb nodes in each interval

% define the function that we wish to approximate
fun = @(x) log(x + 3).*exp(-2*x).*(x.^2).*cos(20*x);

% define the domain over which f is defined
ax = -2; bx = -1;
chebApproxInfo.domain = [ax,bx]; 
chebApproxInfo.ninters = ninters; 
chebApproxInfo.nord = nord; 

% table of cheb weights times function values
wts = mktab_wts(fun, chebApproxInfo);

Nx = 120;
x = linspace(ax,bx,Nx)';
y = chebEval(x,wts,chebApproxInfo);

plot(x,fun(x),'xb',x,y,'or','LineWidth',2);
xlabel x; ylabel y; legend('true','approx');
title(sprintf('Approximation of f(x) = x^2e^{-2x}cos(20x)log(x + 3) via\n order %d Chebyshev polynomials\n on %d subintervals',nord,ninters))
set(gca,'fontsize',15);
end