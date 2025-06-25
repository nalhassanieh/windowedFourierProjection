function [u,ue,tn_sol,n_sol,tSOL_now,sn_sol] = get_solnAndErr(u,ue,tn_sol,n_sol,...
        tSOL_now,x,t,snIdxNow,sn,an,interpNodesIdx_sol,...
        interpAndGLWeights_sol,interpShift_sol,s,M,delta,tol,tSOL,solnType,dataParam,fixSolCnt,n,tgh,sn_sol)
% GET_SOLNANDERR to compute the solution and error (if available)
% 
% INPUTs
%   u,ue: computed and exact solutions
%   tn_sol,n_sol,tSOL_now: solution time grid and index, tSOL_now is a
%   variable to track the time of the solution with respect to the actual
%   time grid
%   x: spatial grid for testing
%   t,tn: current time, and time grid
%   snIdxNow: index of sn at the first time-step
%   sn: the density grid functions
%   an: Fourier Coefficients
%   interpNodesIdx_sol,interpAndGLWeights_sol,interpShift_sol:
%   interpolation parameters needed for the computation of the local part
%   of the solution. 
%   s,M: source location and total num of sources 
%   tol: error tolerance
%   mu,t0: parameters needed for the computation of the manufactured
%   solution for error computation. These variables are held in the
%   dataParam struct for ease of passing in functions and switching between
%   test solutions. 
%   W,delta: width of the window function: W*dt = delta
%   h_temp: paramter needed for free-space projection
%   tSOL: time-step for the test solution grid 
%   addProj: button to get free-space solution 
%   solnType: manufactured or true
%   writerObj: object to save movie 
%   plotMovie: button to plot movie
%   saveMovie: button to save movie
% 
% OUTPUTs
%   u,ue: computed and exact solutions
%   tn_sol,n_sol,tSOL_now: solution time grid and index, tSOL_now is a
%   variable to track the time of the solution with respect to the actual
%   time grid


mu = dataParam.mu; 
t0 = dataParam.t0;

% if(t>=tSOL_now)
if(mod(n-tgh,fixSolCnt) == 0 && t>=0)

    sn_sol(:,(n_sol + 1)) = sn(end,:)';
    
    u(:,(n_sol+1)) = get_uSol(x,t,snIdxNow,sn,an,interpNodesIdx_sol,...
        interpAndGLWeights_sol,interpShift_sol,s,M,delta,tol);

    if(strcmp(solnType,'ms'))
        ue(:,(n_sol+1)) = get_ExactSol(x,t,mu,t0,M,s); 
    else
        ue(:,(n_sol+1)) = get_uIncident(x,t,mu,t0);
    end
    
    tn_sol(n_sol+1) = t;
    n_sol = n_sol + 1;
    tSOL_now = tSOL_now + tSOL;
end

end