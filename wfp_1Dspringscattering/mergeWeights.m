function [wv_out] = mergeWeights(A,w,Wth)
% MERGEWEIGHTS a function to add weights at similar nodes to avoid double
% sums 
% 
% Inputs: 
%   A: some matrix A with some repeated elements 
%   w: weights corresponding to each element of matrix A
% Outputs: 
%   v: the collapsed matrix A with no repeated elements 
%   wv: the sum of weights of repeated elements of A. Each such sum
%   corresponds to the unique value of the element saved in v. 

[m,n] = size(A);
v = unique(A);   % create a vector with unique A elements 
wv = zeros(length(v),1);
for i = 1:m
    for j = 1:n
        idx = find(v == A(i,j));
        wv(idx) = wv(idx) + w(i,j); 
    end
end
 
wv_out = zeros(Wth,1); 
wv_out(v(1):(v(1) + length(v)-1)) = wv; 

end