function x = get_xgrid(ax,bx,src_dmn,Nx,h0,W)
% GET_XGRID a function to get the spatial grid for plotting
% 
% INPUTS: 
%  [ax,bx]: spatial domain
%  [as,bs]: interval where sources are placed
%       Nx: number of grid points in the spatial discretization
%       h0: intial time-step size
%        W: width of the window function
%  addProj: button to turn on free-space projection
%  getSelfConvergence: button to turn on self-convergence study

as = src_dmn(1); bs = src_dmn(2); 

h_temp = h0;

xmin = (ax + 2*h_temp*W); xmax = (bx - 2*h_temp*W);
x = linspace(xmin,xmax,Nx)';
if(xmin>as || xmax<bs)
   xmin,xmax,h0
    error('Sources are outside allowed region; reduce the time-step or narrow down the domain where sources live');
end

end