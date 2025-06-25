function [wts,LebesgueConstant] = get_lmmwts(xpts, npts, nqloc,rlen,dt,delta,winData,computeLebesgueCst)
% GET_LMMWTS gives the weights for the LMM-style treatment of [0,dt]
% integrals or [0,L] integrals where L<dt
%
% wts = get_lmmwts(xpts, npts, nqloc,rlen,phit,blendingTol,dt)
%  performs barycentric interpolation using 'xpts' (a vector of size 'npts')
%  as interpolation nodes, then applies GL quadrature with 'nqloc'
%  quadrature nodes to evaluate the integral over
%  [xpts(npts-1),xpts(npts-1) + rlen].
%  The weights are modulated by the blending function needed for the windowing
%  in the analytic split of the solution. Here 'phi_wts' are cheb weights
%  to compute the blending function.
%   'chebApproxInfo': a struct with fields domain [ax,bx] over which cheb
%     weights are computed; ninters, number of subintervals, and nord,
%     order of polynomial interpolation or number of cheb nodes.
%  'computeLebesgueConstant' is a button to compute Lebesgue constant if
%  needed (set to one).
%
% Note: function rewritten from Fortran function provided by Leslie
%  Greengard

if nargin==0, test_getlmmwts; return; end

rintmatloc = zeros(nqloc, npts); % Allocate memory for integration matrix

% Obtain Gauss-Legendre nodes and weights -> on [-1,1] interval
gl = glwt_prep(nqloc);
glnodest = gl.x0; 
glweightst = gl.w0; 

% precompute the barycentric weights
X=repmat(xpts,1,npts);
w = 1./prod(-X+X.'+eye(npts),1);

% Loop over intervals (xpts(1),..,xpts(npts))
aa = xpts(npts-1); % here the interval of integration [aa,bb] = [xpts(npts-1),xpt(npts)]

% GL nodes adjusted to interval
xt = aa + rlen*(glnodest + 1.0)/2.0;
phi = generalwindow((dt - xt)/delta,winData);
phimult = (1 - phi);

glwloc = rlen*glweightst.*phimult;

tol = 1.0e-20;
for iquad = 1:nqloc % loop over each GL node
    rr = prod((xt(iquad) - xpts));
    rintmatloc(iquad, :) = rr*w./(xt(iquad) - xpts');
    rintmatloc(iquad,(abs(xt(iquad) - xpts) < tol)) = 1;
end

% Compute the Lebesgue constant
if(computeLebesgueCst == 1)
    LebesgueConstant = norm(rintmatloc,inf);
else
    LebesgueConstant = NaN;
end

% Compute weights: GL quadrature applied to each lj
wts = glwloc'*rintmatloc; 

end % end of get_lmmwts function

function test_getlmmwts
clf;
ax = 0; bx = 1;
h = 1e-4;
xpts = ax:h:bx;
rlen = h;
tol = 1e-12;    % error tolerance
gam = .5;       % pad gam*Nyquist for window
theta = log(1/tol);
W = ceil(2*theta/(pi*gam));
[phit1,phitt1] = window(theta,W);
chebApproxInfo = struct('ninters',12,'nord',5,'domain',[-36,100]);
[phi_wts,~,~]= prepWindowChebWts(phit1,phitt1,chebApproxInfo,tol);

% npts = 7;
% nqloc = npts;
% [wts,LebesgueConstant] = get_lmmwts(xpts, npts, nqloc,rlen,phi_wts,h,chebApproxInfo);
% bar(wts);
% title(sprintf('Order %2d LMM weights for grid points\n with separation %1.2e',npts,h));

p = 8;
for i = 1:p
    npts = 2*i;
    nqloc = npts;
    [wts,LebesgueConstant] = get_lmmwts(xpts, npts, nqloc,rlen,phi_wts,h,chebApproxInfo);

    subplot((p/2),2,i)
    bar(wts);
    str = '\Lambda';
    title(sprintf('order $%2d$, $%s_{%d} = %1.1f$',npts,str,npts,LebesgueConstant),'Interpreter','latex');
end

sgtitle(sprintf('LMM weights for grid points\n with separation %1.2e',h))
saveas(gcf,'/Users/nalhassanieh/Desktop/stabilityResults/LMM_wts2','epsc');
end
