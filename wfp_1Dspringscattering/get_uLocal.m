function ul = get_uLocal(x,t,snIdxNow,sn,interpNodesIdx_sol,interpAndGLWeights_sol,interpShift_sol,s,M,delta)
% GET_ULOCAL computes the local part of the solution 
%
% ul = get_uLocal(x,t,snIdxNow,sn,interpNodesIdx_sol,...
%                 interpAndGLWeights_sol,interpShift_sol,s,M,delta)
% returns the local part of the solution at each value of the vector 'x' at
% time 't'. The function takes the unifrom time-grid 'tn', 'delta' ([0,delta]
% is the compact support of the window), 'M' source locations in the vector
% s = [s(1),...,s(M)], the computed density values sn.
% snIdxNow: index of the density at the current time step. 
% interpNodesIdx_sol: interpolation nodes used to evaluate the density at
%  GL nodes
% interpAndGLWeights_sol: interpolation/integration weights needed to
%  evaluate the local integral
% interpShift_sol: initial shift of interpolation nodes (to undo the shift)

% Pick interpolation points between [startTime,endTime], where startTime
% and endTime lie on the time-grid
Nx = length(x); 
    ul = zeros(Nx,1);
    for i = 1:Nx % for each x value
        for j = 1:M % for each source
            L = abs(x(i) - s(j)); 
            if(L<delta && (t - L)>0)
                snj = sn(:,j); 
                % perform the local integration using GL and barycentric
                % interpolation
                 shift = (snIdxNow - 1) + interpShift_sol; 

                 shiftedInterpNodesIdx = interpNodesIdx_sol{i,j} + shift;
                 snIdx = shiftedInterpNodesIdx + 1; %snIdx(snIdx<=0) = 1;
                 snvec = sum(interpAndGLWeights_sol{i,j}.*snj(snIdx),2);

                 I = sum(snvec);
                 ul(i) = ul(i) + I;
            end
        end
        ul(i) = 0.5*ul(i);
    end

end

% Note replaced for loops by matvecs. Original:
% snvec = zeros(P,1);
% for k = 1:P  % for each GL node calculate interpolated value at GL node
%     shiftedInterpNodesIdx = interpNodesIdx_sol{i,j}(k,:) + shift;
%     snIdx = shiftedInterpNodesIdx + 1; snIdx(snIdx<=0) = 1;
%     snvec(k) = interpAndGLWeights_sol{i,j}(k,:)*(sn{j}(snIdx));
% end
% I = sum(snvec); 
