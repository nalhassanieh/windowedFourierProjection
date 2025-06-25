function g = get_gvalueFromIncomingField(t,s,j,beta,mu,t0)
% GET_GVALUEFROMINCOMINGFIELD get data from incoming field
%
% g = get_gvalueFromIncomingField(t,s,j,beta) returns the data g for each spring
%  scatterer condition at a time 't'. 's' is a vector containing the source
%  locations s = [s(1),...,s(M)], where 'M' is the number of sources. 
%  'beta' is a factor of the spring constant, and 'j' represents the index 
%  of the current spring scatterer. 'mu' and 't0' are parameters needed for
%  the incident wave. 

x = s(j); 
uin = get_uIncident(x,t,mu,t0);
g = beta*uin; 

% testing 
% g = beta*cos(10*(s(j) - t - t0)); 

end