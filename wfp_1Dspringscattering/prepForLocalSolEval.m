function [interpNodesIdx_sol,interpAndGLWeights_sol,interpShift_sol] = prepForLocalSolEval(x,dt,gl,s,M,P,W,order,winData)
% PREPFORLOCALSOLEVAL precomputes parameters for local solution evaluation 
%
% [interpNodesIdx_sol,interpAndGLWeights_sol,interpShift_sol] = ...
%  prepForLocalSolEval(x,dt,delta,gl,s,M,phi_wts,P,W,order,chebApproxInfo)
%      prepares interpolation nodes and weights needed for the evaluation
%      of the local part of the solution. The function takes 'x' the
%      spatial grid, 'dt' the time-step, 'delta = W*dt' the width of the window
%      function 'gl' GL nodes and weights on [-1,1], 's' a vector
%      carrying the location of the sources, 'M' the number of sources, the
%      blending function cheb weights 'phi_wts' and the number of GL nodes 'P'.
% 
% NOTE:'chebApproxInfo': a struct with fields domain [ax,bx] over which 
%      cheb weights are computed; ninters, number of subintervals, and nord, 
%      order of polynomial interpolation or number of cheb nodes     
%
% OUTPUT: 
%   interpNodesIdx_sol: indices of interpolation nodes
%   interpAndGLWeights_sol: weights needed to perform the integration of
%     the local part of the soluiton 
%   intepShift_sol: how many grid points are the interpolation points
%     shifted. We undo the shift when we use these nodes in the
%     time-stepping loop. 

delta = W*dt; 

Nx = length(x); % number of points in the spatial grid

% allocate space for variables
interpNodesIdx_sol = cell(Nx,M); interpNodesIdx_sol(:) = {zeros(P,order)}; 
interpAndGLWeights_sol = cell(Nx,M); interpAndGLWeights_sol(:) = {zeros(P,order)};

interpShift_sol = -(W-1); % shift index of the interpolation nodes (nodes shifted forward by delta)

for i = 1:Nx % for each x vlue
    for j = 1:M % for each source
        L = abs(x(i) - s(j));
        if(L<delta)

            % prepare the GL nodes eta needed for window factor 
            [eta,~] = glwt(L,delta,gl);

            % Get the blending function at eta
            phiL_to_delta_sol = generalwindow((eta./delta),winData);

            % actual interval of integration 
            % t1 = t - delta; 
            % t2 = t - L;
             
            % interval of integration (will be shifted later)
            t1 = 0;
            t2 = delta - L;

            % Get the GL nodes and weights in the interval of integration 
            [glnodes,glwts] = glwt(t1,t2,gl);

            % update the glwts by the window factor
            glwts = glwts.*(1 - phiL_to_delta_sol(end:-1:1));

            % perform Barycentric interpolation at the GL nodes and get
            % weights 
            interpWeights_sol = zeros(P,order); 
            endTime = t2; 
            for k = 1:P
                [interpNodes,interpNodesIdx_sol{i,j}(k,:)] = get_InterpNodes(glnodes(k), dt, order,endTime);
                interpWeights_sol(k,:) = barycentricInterp(interpNodes, order, glnodes(k));
            end

            % Multiply the interpolation weights by the GL weights 
            interpAndGLWeights_sol{i,j} = interpWeights_sol.*glwts;
        end
    end
end

end