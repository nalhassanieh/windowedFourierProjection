function L = antiderivmat_1d(s)
% ANTIDERIVMAT_1D   matrix from nodes in 1D to antiderivatives at nodes
%
% L = antiderivmat_1d(s) returns (N-1)*N matrix taking values on s,
%  a list of N nodes, to their integrals from the first node to each of the
%  other N-1 nodes in turn. Ie, the antiderivative with constant chosen so that
%  its value at the first node would be zero. Nodes must be spaced in a sensible
%  way. length(s) should not exceed around 30 for stability reasons.
%
% Notes: 1) should match Leslie's get_integratemat().
%   2)  Computed in Helsing style, with centering and scaling for stability.

% Barnett 8/13/24.

if nargin==0, test_antiderivmat_1d; return; end

cen = (max(s)+min(s))/2; hwid = (max(s)-min(s))/2;   % affine map s to [-1,1]
s = (s-cen)/hwid;
n = numel(s); s = s(:);       % col vec
V = ones(n); for j=2:n, V(:,j) = V(:,j-1).*s; end   % Vandermonde (polyval) mat
U = diag(s)*V*diag(1./(1:n));          % mat evaluating an antideriv of poly
L = (V'\U')';    % backwards-stable way to solve for it (Helsing)
L = L(2:end,:) - L(1,:);      % adjust const to zero at first node
L = L*hwid;                   % unscale

%%%%%%
function test_antiderivmat_1d
off = 4.3; sc = 1.7;          % test centering and scaling
x = sc*linspace(-1,1,16)' + off;  % nodes
%x = sc*gauss(16) + off;       % nodes
f = @(x) sin(0.8*x + 0.7);    % the antiderivative
fp = @(x) 0.8*cos(0.8*x + 0.7);   % the input func
Fex = f(x(2:end))-f(x(1));    % exact ans at nodes 2...N, col vec
L = antiderivmat_1d(x);
F = L * fp(x);                % hit L against vec of func values
fprintf('max abs err for antideriv on nodes : %.3g\n',max(abs(F - Fex)))
fprintf('[mat inf-norm = %.3g;   max element size = %.3g]\n',norm(L,inf),max(abs(L(:))))
