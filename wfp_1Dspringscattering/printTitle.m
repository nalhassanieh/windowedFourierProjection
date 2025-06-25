function printTitle(tFinal, tol, W, P, M, typNumOfNeighbors, order, solnType)
% PRINTTITLE  print a title for the driver file indicating important
% problem parameters
% 
% printTitle(tFinal, tol, W, P, gam, order, src)
%  Displays important problem paramters in the command window
% 
% Input values: 
%  tFinal: final time
%  tol: tolerance for the window function
%  W: W*dt is the support of window function
%  P: number of GL nodes usually set to W
%  gam: pad gam*Nyquist for window
%  order: order of accuracy 
%  M: number of spring scatterers 
%  solnType: type of testing solution 

fprintf('wave: t = %3.0f, tol = %1.0e, W = %d, P = %d, M = %d, typNumNeighbors = %d, order = %d, solnType = %s\n',...
    tFinal, tol, W, P, M, typNumOfNeighbors, order, solnType); 

end