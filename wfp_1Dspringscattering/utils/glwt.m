function [x,w] = glwt(a,b,gl)
%GLWT convert GL nodes and weights over [-1,1] to ones over [a,b]
%
%[x,w] = glwt(a,b,x0,w0) returns GL nodes and weights x and w over [a,b]
% given GL nodes and weights x0,w0 over [-1,-1]

% get nodes and weights
x0 = gl.x0; 
w0 = gl.w0; 

% Linear map from[-1,1] to [a,b]
x=(a*(1-x0)+b*(1+x0))/2;      

% Compute the weights
w=(b-a)*w0; 

end