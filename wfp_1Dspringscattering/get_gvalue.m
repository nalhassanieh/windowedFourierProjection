function g = get_gvalue(t,s,j,sig,beta,M,mu,t0,addProj)
% GET_GVALUE get the data at a given time t
%
% g = get_gvalue(t,s,j,sig,beta,M,mu,t0) returns the data g for each spring
%  scatterer condition at a time 't'. 's' is a vector containing the source
%  locations s = [s(1),...,s(M)], where 'M' is the number of sources. 'sig'
%  is a cell that holds manufactured density values (Gaussians) where
%  sig{j} = @(t) exp(mu(j)*(t - t0(j))). 'beta' is a factor of the spring
%  constant, and 'j' represents the index of the current spring scatterer

xval = s(j); % location at the current source
U = get_ExactSol(xval,t,mu,t0,M,s,addProj);  % get the manufactured solution at xval
g = -sig(t,j) - beta*U;

end