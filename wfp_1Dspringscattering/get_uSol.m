function u = get_uSol(x,t,snIdxNow,sn,an,interpNodesIdx_sol,interpAndGLWeights_sol,interpShift_sol,s,M,delta,tol)
% GET_USOL gets the solution by forming local and history parts
% 
% u = get_uSol(x,t,snIdxNow,sn,an,K,interpNodesIdx_sol,...
%              interpAndGLWeights_sol,interpShift_sol,s,M,delta) 
%  returns the solution u by evaluating the local and the history parts.
%  The function requires the spatial grid 'x', the current time 't', the
%  unifrom time grid 'tn', the Fourier coefficients 'an', the frequency
%  grid 'K', 'delta' the width of the window function used in the
%  local/history split, the location of sources 's' s = [s(1),...,s(M)] 
%  where 'M' is the number of sources, 'sn' is the computed density.
% 
% Other input values: 
% snIdxNow: index of the density at the current time step. 
% interpNodesIdx_sol: interpolation nodes used to evaluate the density at
%  GL nodes
% interpAndGLWeights_sol: interpolation/integration weights needed to
%  evaluate the local integral
% interpShift_sol: initial shift of interpolation nodes (to undo the shift)

% get the local part of the solution 
ul = get_uLocal(x,t,snIdxNow,sn,interpNodesIdx_sol,interpAndGLWeights_sol,interpShift_sol,s,M,delta);

% get the history part of the solution 
uh = get_uHistory(an,x,tol);

% compute the solution 
u = uh + ul;

end