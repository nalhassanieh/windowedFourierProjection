function TEXP = chexfcdir(F, N)
% CHEXFCDIR computes coefficients of Chebyshev expansion
% of function tabulated at classical Chebyshev nodes.
%
% INPUT:
% F    = function values at Chebyshev nodes in increasing order
%        on the interval (-1,1)
% N    = number of Chebyshev nodes
%
% OUTPUT:
% TEXP = array of expansion coefficients of (N-1)st degree
%        interpolant
% 
% NOTE: code adjusted from code provided by Leslie Greengard. 

if nargin == 0, test_chexfcdir; return; end

F = F(:); % ensure the input is a column vector

FAC = 2/N;       % Factor to normalize
TEXP = zeros(N, 1); % Initialize TEXP array

J = (1:N)'; % cheb node index
for I = 1:N % loop for each cheb coefficient
    TEXP(I) = sum(FAC * F(J).*cos((J - 0.5)*(I - 1)*pi./N));
    TEXP(I) = -TEXP(I)*(-1)^I; % Apply the sign change
end
TEXP(1) = TEXP(1) / 2;  % Adjust the first coefficient

end

function test_chexfcdir
N = 16;
[x, ~, ~,  ~,  ~,  ~] = chnodc(-1, 1, N);

F = cos(x - 3).*exp(-2*x)';

TEXP = chexfcdir(F, N);
disp('Array of expansion coefficients of (N-1)st degree interpolant');
disp(TEXP);
end