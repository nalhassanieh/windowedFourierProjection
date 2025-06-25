function wts = mktab_wts(fun, chebApproxInfo)
% MKTAB_WTS make a table of function values over cheb nodes in subintervals
% of domain [ax,bx]
% 
% INPUTs: 
%   fun: function handle
%   chebApproxInfo.domain: domain [ax,bx] overwhich wts are computed. 
%   chebApproxInfo.ninters: number of subintervals
%   chebApproxInfo.nord: order of polynomial interpolation; number of cheb 
%         nodes
%
% OUTPUT
%   wts(norder,ninters): table of cheb weights in each of
%   the subintervals where cheb poly approximation is applied
%
% NOTE: code adjusted from code provided by Leslie Greengard. 

if(nargin == 0), test_mktab_wts; return; end

nord = chebApproxInfo.nord; 
ninters = chebApproxInfo.ninters; 
domain = chebApproxInfo.domain; 
ax = domain(1); bx = domain(2);  

% Define function
% fun = @(x) ...;

% Preallocate the output matrix
wts = zeros(nord, ninters);

h = (bx - ax)/ninters; 

for iint = 1:ninters
    A = ax + (iint - 1)*h;
    B = ax + iint*h;

    % Get Chebyshev nodes and sine values
    [CHPTS,~,~,~,~,~] = chnodc(A, B, nord);
    wts(:, iint) = chexfcdir(fun(CHPTS), nord);
end
end

function test_mktab_wts
ninters = 10;  % Number of intervals
nord = 5;      % Order for Chebyshev nodes
fun = @(x) cos(x - 3).*exp(-2*x);

ax = -1; bx = 1;
chebApproxInfo.nord = nord; 
chebApproxInfo.ninters = ninters; 
chebApproxInfo.domain = [ax,bx]; 

wtsTab = mktab_wts(fun, chebApproxInfo);

disp('Function at cheb nodes in subintervals:');
disp(wtsTab);
end