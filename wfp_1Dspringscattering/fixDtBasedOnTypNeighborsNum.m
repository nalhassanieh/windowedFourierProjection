function [dt,typNumOfNeighbors] = fixDtBasedOnTypNeighborsNum(maxNumNeighbors,src_dmn,W,M,dt)
% FIXDTBASEDONTYPNEIGHBORSNUM adjusts dt to make sure num of local
%  neighbors is less than a maximum value
% 
% INPUTS:
% maxNumNeighbors: maximum number of neighbors
% [as,bs]: domain in which the sources are located 
% W: W*dt is the compact support of the window function
% M: number of spring scatterers
% dt: time step (first estimate) 
% N: 2pi/dt: number of Fourier modes
% tFinal: final time. 
%
% OUTPUTS:
% dt: adjusted timestep
% N: adjusted number of Fourier modes
% typNumOfNeighbors: typical number of neighbors expected on average for
%  each source.

as = src_dmn(1); bs = src_dmn(2); 

if(isnan(maxNumNeighbors))
    maxNumNeighbors = ceil(sqrt(M*log(2*W*M))); % choose to balance otherContrib and finufft
end

typNumOfNeighbors = ceil(2*W*dt*M/(bs - as));

if(typNumOfNeighbors>maxNumNeighbors)
    dt = (bs - as)*maxNumNeighbors/(2*W*M);
    fprintf('time step changed to dt = %1.2e to ensure maximum number of neighbors %d\n',dt,maxNumNeighbors);
    typNumOfNeighbors = maxNumNeighbors; 
end


end