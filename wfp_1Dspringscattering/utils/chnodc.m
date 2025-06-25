function [CHPTS, SINCH, U, V, U1, V1] = chnodc(A, B, M)
% CHNODC Constructs Chebyshev nodes and mapping coefficients.
%
%   [CHPTS, SINCH, U, V, U1, V1] = CHNODC(A, B, M)
%
%   This function calculates classical Chebyshev nodes for the
%   interval [A, B] and the coefficients for linear mappings
%   between the intervals [-1, 1] and [A, B].
%
%   INPUT PARAMETERS:
%   A - Lower bound of the interval.
%   B - Upper bound of the interval.
%   M - Number of Chebyshev nodes to generate.
%
%   OUTPUT PARAMETERS:
%   CHPTS - Array of Chebyshev nodes in the interval [A, B].
%   SINCH - Array containing U * sin(theta) for the Chebyshev points.
%   U - Coefficient for mapping from [-1, 1] to [A, B].
%   V - Coefficient for mapping from [-1, 1] to [A, B].
%   U1 - Coefficient for mapping from [A, B] to [-1, 1].
%   V1 - Coefficient for mapping from [A, B] to [-1, 1].
%
%   The kth Chebyshev node is computed using:
%       CHPTS(k) = U * cos((2k - 1) * pi / (2M)) + V
%
%   Example:
%       A = -1;
%       B = 1;
%       M = 5; % Number of Chebyshev nodes
%       [CHPTS, SINCH, U, V, U1, V1] = chnodc(A, B, M);
%
% NOTE: code adjusted from code provided by Leslie Greengard. 

if(nargin == 0), test_chnodc; return; end

% Construct the scaling parameters
U = (B - A) / 2;
V = (B + A) / 2;
U1 = 2 / (B - A);
V1 = 1 - U1 * B;

% Preallocate the arrays for Chebyshev nodes and sin values
CHPTS = zeros(1, M);
SINCH = zeros(1, M);

% Construct the Chebyshev nodes and corresponding SIN array
K = 1:M;
CHPTS(M - K + 1) = U * cos((2 * K - 1) * pi / (2 * M)) + V;
SINCH(M - K + 1) = U * sin((2 * K - 1) * pi / (2 * M));

end

function test_chnodc
clf;

A = -2;
B = 1;
M = 20; % Number of Chebyshev nodes
[CHPTS, SINCH, U, V, U1, V1] = chnodc(A, B, M);

plot(CHPTS,0,'*r','markerSize',10);
end