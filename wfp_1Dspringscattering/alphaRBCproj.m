function [an,bn,info] = alphaRBCproj(an,bn,phi,phip)
% ALPHARBCPROJ  apply projector to Fourier coeff rep for 1D free space BCs
%
% [an,bn] = alphaRBCproj(an,bn,phi,phip) applies RBCs to the alpha_n and beta_n
%  vectors (set of Fourier coefficients of u and v=u_t, respectively).
%  phi is a length-w vector giving a (roll-up) window evaluated on
%  the Nyquist grid, and phip the evaluation of its derivative (scaled as if
%  the grid were unit sized). w sets the tolerance, and phi and phip could
%  be as in the history split. The cost is O(N log N), where
%  N=length(an)=length(bn) is the number of coeffs (N/2 the max frequency).
%
% Usage:
%  The user's computational domain is Omega = (-pi+2wh,pi-2wh), where
%  h = 2pi/N is the Nyquist grid spacing implied by the number of coefficients,
%  for the periodic domain [-pi,pi). The user simulation must live in Omega.
%  [-pi,-pi+2wh] and [pi-2wh,pi) are padding
%  regions needed to apply the RBCs, and waves in them will be destroyed.
%  Thus the user sacrifices 4w of the N points to the RBCs, in order that
%  the Fourier series well approximates the free-space solution in Omega.
%  To maintain outgoing Fourier coefficients, this RBC projection must be
%  reapplied in tRBC time units, for wh <= tRBC <= 2wh.
%  The lower limit is needed so that traveling waves have cleared the
%  windowing region (otherwise the window hits them multiple times and
%  aliasing error occurs). The upper limit is so that waves do not wrap around
%  into the other transition region, since the assumption is made that the
%  waves in transition region (eg (pi-2wh,pi-wh)) are out-going.
%  Even without this assumption, spectral differentiation would be needed
%  to recover v in the padding regions, but waves would be hit >1 times by
%  the window, causing aliasing error. (For t > 3wh,
%  plain old periodic pollution would also occur back in Omega).
%  The phi and phip should be as in the history split (as given by window.jl).
%
% [an,bn,info] = alphaRBCproj(an,bn,phi,phip) also prints & returns debug info.
%
% When called with no arguments, a self-test is done.
%
% The method is:
%  1) evaluate u and v on the h grid, via the FFT from the coeffs an and bn.
%  2) multiplication of u by phi-windows in the padding regions (which
%     broadens the spectral bandwidth, hence gam is needed),
%  3) building v=u_t=+-u_x on the grid points in the padded regions to give
%     outgoing 1D wave equation solutions there.
%  4) inverse FFT to the output coefficients an, bn.
%
% Notes:
%  i) N even for now.
%  ii) The total padding region could easily be shrunk from 4w to 3w points.
%     The former gives a bit more flexibility in the projection interval.
%  iii) check Fourier series defn, sign, prefactor, may not matter.

% Barnett 9/8/24
if nargin==0, test_alphaRBCproj; return; end
verb = (nargout>2);         % verbosity
N = length(an); assert(length(bn)==N)
w = length(phi); assert(length(phip)==w)
assert(4*w<N)
h = 2*pi/N;

% 1)
u0 = fft(an(:)); v0 = fft(bn(:));  % x on [0,2pi), so padding is in middle

% 2)
j = N/2-2*w+(1:4*w)';   % indices to bleach out (combine R and L padding)
phi=phi(:); phip=phip(:);   % force col vecs
phij = [flipud(phi); zeros(2*w,1); phi];   % R then L padding
pou = ones(N,1); pou(j) = phij; info.pou=pou;  % save the POU, full [0,2pi) grid
u = u0 .* pou;
v = v0 .* pou;     % the naive step 3, alone would cause O(1/w) reflection.
% This assumes eg, waves in (pi-2wh,pi-wh) are R-going

% 3) apply prod rule, overwrite v = (phi.u0)_x = phi.u0_x + phi'.u0 in padding,
% noting that phi.u0_x is already accounted for in v.
phipj = (1/h) * [flipud(phip); zeros(2*w,1); phip];   % takes d/dx or -d/dx
v(j) = v(j) + phipj.*u0(j);                 % add term +- phi'.u0 in prod rule

% 4)
an = reshape(ifft(u),size(an)); bn = reshape(ifft(v),size(bn));


%%%%%% helpers for test
function [u,v,xg] = show_uv(an,bn)  % plot Fourier rep in real space on Nyq grid
N = length(an);
xg = -pi + (0:N-1)/N*2*pi;
u = fftshift(fft(an)); v = fftshift(fft(bn));   % shift since grid starts -pi
plot(xg, real([u;v]), '.-')
a=axis; a(1:2)=[-pi,pi]; axis(a); drawnow

function [an, bn] = propexact(t, kg, an0, bn0)  % evol by t the coeffs on k grid
ak=abs(kg); s=sin(ak*t); c=cos(ak*t);         % (taken from tryEFrep.jl)
sok=s./ak; sok(kg==0) = t;                    % note all act on vectors
an =      c.*an0 + sok.*bn0;
bn = -s.*ak.*an0   + c.*bn0;

%%%%%%
function test_alphaRBCproj   % test can evolve a source-free WE in free space
addpath ./utils
verb = 1;         % 1,2,.. for debug tests
tol = 1e-12;
gam = 0.5;          % fractional exceeding of Nyquist of the sigma

theta = log(1/tol);        % set up window on regular grid in Nour's way
w = ceil(2*theta/(pi*gam));
phipfun = window(1.0,theta,w);   % hack dt=1  (not so happy window needs dt)
g = (1:w) - 0.5;           % not sure of off-by-one shifts here
phip = phipfun(g)';    % col vec
phi = blending(phipfun, g, tol);
if verb>1
    figure; plot(g, [phi, phip], '-+'); axis tight; xlabel('gridpoints');
    legend('\phi','\phi'''); drawnow
    % check that phip is really deriv of phi on the same grid...
    foldphipe = [phip;flipud(-phip)];  % append flipped so can do periodic diff
    foldphip = (pi/w) * perispecdiff([phi;flipud(phi)]); % 2pi/2w rescales t
    fprintf('max err in phip on grid: %.3g\n', norm(foldphip-foldphipe,inf))
end

N = 300;  assert(mod(N,2)==0)     % grid size N even
x0 = 0.8;               % test initial condition center loc
k0 = N/sqrt(8*theta);      % k-width of Gaussian which hits tol by N/2
kg = [0:N/2-1, -N/2:-1];  % freq k grid, fft ordering (may differ from Nour)
c0 = 2/(sqrt(2*pi)*k0);    % Gaussian height 2 (so splits into packets of 1)
an = c0 * exp(-0.5*(kg/k0).^2) .* exp(1i*x0*kg);  % bump at x0
bn = 0*an;
h = 2*pi/N;    % implied grid (used inside proj only, plus for plotting)
tRBC = 1.0*w*h;     % may check fails either side of allowed range [1,2] :)
fprintf('%.3g of N is in Omega; tRBC=%.3g\n',1-4*w/N,tRBC)
T = 5*pi;
nt = ceil(T/tRBC)+1; tg=tRBC*(0:nt-1); uxt = nan(N,nt);
info.pou = nan(1,N);
if verb, figure; end
for n=1:nt, t = tg(n);        % t step
    if verb, subplot(2,1,1); [u,v,xg] = show_uv(an,bn);
        if verb>1, hold on; plot(xg,fftshift(info.pou),'k-'); hold off; end
        title(sprintf('u and v: t=%6.2f',t));
        vline([-1 1]*(pi-2*w*h)); vline([-1 1]*(pi-w*h));  % POU on/off pts
        uxt(:,n) = u(:);
        subplot(2,1,2); imagesc(xg,tg,log10(abs(uxt))'); axis xy;
        caxis([-12 0]); xlabel x; ylabel t; colorbar; title('log_{10}|u(x,t)|')
        vline([-1 1]*(pi-2*w*h)); vline([-1 1]*(pi-w*h));  % POU on/off pts
        drawnow;
    end
    [an,bn,info] = alphaRBCproj(an,bn,phi,phip);  % comment to switch off
    [an, bn] = propexact(tRBC, kg, an, bn);   % we plot *after* prop
    pause(.5);  % anim already too slow anyway :(
end



