function [phit,phitt] = window(b,w)
% WINDOW1 generates the Kaiser Bessel window function and its time
% derivative
%
% [phit,phitt] = window(b,w)
%  returns the Kaiser Bessel window function phit and its time derivative
%  phitt whose support is w. 
%
% Input values: 
%  b = log(1/tol), where tol phit(0) = phit(delta) = tol
%  w : the width of the support of phit
%
% Note: fuNction rewritten from a Julia code provided by Alex Barnett

cen = w/2.0; % center of the window function
prefac = (1/w) * b/(sinh(b)); 
phit = @(x) prefac*besseli(0,b*sqrt(1-((2/w)*(x-cen)).^2));
phitt = @(x) prefac*besseli(1,b*sqrt(1-((2/w)*(x-cen)).^2)).*(-4*b*(x - cen))./((w^2)*sqrt(1-((2/w)*(x-cen)).^2));

end