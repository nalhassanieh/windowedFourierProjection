function U = get_ExactSol(x,t,mu,t0,M,s)
% GETEXACTSOL gives the exact solutions for a Gaussian densities
%
% U = get_ExactSol(x,t,mu,t0,M,s,addProj)
%  Compute an exact solution for testing by choosing Gaussian densities
%  with at locations 'x' (may be an input vector) and time 't' (one fixed 
%  time value). The Gaussians take the form
%  sig(t) = exp(-mu*(t - t0)^2) where 'mu' and 't0' are vectors of size 'M'
%  where 'M' is the number of sources and 's' represents the source
%  locations s = [s(1),...,s(M)]
%
% If addProj == 1, the exact solution computed is the one corresponding to
%  the solution of the free-space problem
% If addProj == 0, the exact solution is the one corresponding to the
%  solution of the periodic problem

if nargin == 0, test_getExactSol; return; end

Ng = length(x);
U = zeros(Ng,1);

sqrtpi = sqrt(pi);
sqrtmu = sqrt(mu);

%%% Without the images
for i = 1:Ng % for each x value
    for m = 1:M % for each source
        if(t>abs(x(i) - s(m))) % check needed for integration
            A = t - abs(x(i) - s(m));
            U(i) = U(i) + (sqrtpi/(4*sqrtmu(m)))*(erf(sqrtmu(m)*t0(m)) - erf(sqrtmu(m)*(t0(m) - A)));
        end
    end
end

end % end function

function test_getExactSol
ax = -pi; bx =pi; Nx = 50;
as = -pi; bs = pi;
dt = 1e-1;
t = 0:dt:4*pi; Nt = length(t);
M = 100; Ns = M;
x = linspace(ax,bx,Nx);
s = linspace(as,bs,Ns); s = s(randperm(Ns,M));
mu = linspace(10,50,M);
t0 = linspace(0.5,t(end)-0.5,M);

u = zeros(Nx,Nt);
for n = 1:Nt
     u(:,n) = get_ExactSol(x(:),t(n),mu,t0,M,s);
end

u = real(u);

imagesc(x,t,log10(abs(u))'); axis xy;
% hold on; plot(s,0,'|r','MarkerSize',10); hold off;
a=axis; a(1:2)=[-pi,pi]; axis(a);
xlabel x; ylabel t; colorbar; title('u(x,t)');

% plot(x,U,'LineWidth',2);
% title(sprintf('solution evaluated at t = %1.2f',t));
% xlabel('x'); ylabel('u'); axis tight;

end % end the test function