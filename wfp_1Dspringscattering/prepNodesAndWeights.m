function [implicitMat,A_sp,deltaInteractions,dtInteractions] = ...
    prepNodesAndWeights(tn,dt,beta,gl,s,src_dmn,snIdxNow,P,M,W,order,...
    timeLevels,typNumOfNeighbors,winData)
% PREPNODESANDWEIGHTS to prepare nodes and weights for the numerical
% evaluation of integrals
% 
% INPUTS:
%  tn,dt: time grid (array) with time step dt
%  phi_wts: array of weights to evaluate the blending function using
%   Chebyshev polynomial
%  beta: an array containing spring constants
%  gl: GL nodes gl.x0 and weights gl.w0 over the interval [-1,1]; 
%  s: a vector holding the locations of the spring sources
%  snIdxNow: index of density grid function at the current step. For
%   example, the density functions are arranged in the following order
%   -tgh, ..., -1,0,1, where tgh corresponds to temporal ghost points,
%   snIdxNow corresponds to the index of 0. 
%  P: number of GL nodes 
%  M: number of sources
%  W: width of the window function (i.e. support of the window function up
%   to tol)
%  chebApproxInfo: information needed for the Chebyshev polynomial
%   approximation of window and blending functions
%  timeLevels: the number of time levels stored for the density functions;
%   this is equal to W + order to account for centered barycentric
%   interpolation at the left end of the time axis
%  typNumOfNeighbors: estimated average typical number of neighbors for
%   spring sources 
%
% OUTPUTS:
%  implicitMat: implicit sparse matrix to solve for the density values at
%   the new time
%  A_sp: a sparse matrix to evaluate integrals numerically using a
%  combination of barycentric interpolation and Gauss-Legendre quadrature

% domain of sources
as = src_dmn(1); bs = src_dmn(2); 

% support of the window function
delta = W*dt;

% prepare the implicit matrix
implicitMat = speye(M);

% get the LMM weights for the treatment of selfLocal1
selfLocal1_wts = get_lmmwts(tn(1:order), order, order,dt,dt,delta,winData,0);
coef = -(1 + (selfLocal1_wts(end).*(beta./2))); % coefficient of the sigma at (n+1)
implicitMat = coef.*implicitMat;

% GL nodes and weights for selfLocal2 and blending function on [dt,delta]
[eta,~] = glwt(dt,delta,gl);
phidt_to_delta = generalwindow((eta./delta),winData);

% GL nodes and weights for the selfLocal2 implementation and
% interpolation weights
[selfLocal2_glnodes,selfLocal2_glwts] = glwt(0,(delta-dt),gl);
selfLocal2_interpWeights_temp = zeros(P,order);
interpNodesIdx_temp = zeros(P,order);
for kk = 1:P % for each GL node calculate interpolated value at GL node
    [interpNodes,interpNodesIdx_temp(kk,:)] = get_InterpNodes(selfLocal2_glnodes(kk), dt, order,(delta-dt));
    selfLocal2_interpWeights_temp(kk,:) = barycentricInterp(interpNodes, order, selfLocal2_glnodes(kk));
end
% multiply the weights by GL weights and window factor
selfLocal2_interpWeights_temp = selfLocal2_interpWeights_temp.*(1 - phidt_to_delta(end:-1:1)).*selfLocal2_glwts;

% shift the nodes to fit the position of each node
interpNodesIdx_temp = interpNodesIdx_temp+ snIdxNow + 1 - W;

% merge the weights of similar interpolation nodes
[selfLocal2_interpWeights] ...
    = mergeWeights(interpNodesIdx_temp,selfLocal2_interpWeights_temp,timeLevels);

% create bins to avoid checks over the full set of sources
% each box is delta big, so check sources in the current box and other
% surrounding boxes
nboxes = ceil((bs - as)/delta); % number of boxes
snew = (s - as)./(bs - as);     % transform s so it is on [0,1] for assign function
[ioffst, ibox, isradr,icnt] = assign(nboxes, snew, M); % assign sources to boxes

% Fill in a sparse matrix with weights to evaluate within the time loop
% The product of the sparse matrix with the density vector evaluates
% numerical integrals needed in the solver
nz = M*timeLevels*typNumOfNeighbors;
i_sp = ones(1,nz);
j_sp = ones(1,nz);
a_sp = zeros(1,nz);

cnt = 1;
dtInteractions = zeros(M,1); % set up counts for delta interactions for testing
deltaInteractions = zeros(M,1); % set up counts for dt interactions for testing
for j = 1:M
    adr = ibox(j);
    sourcesNearby_left = []; sourcesNearby_right = [];
    sourcesNearby_center = isradr(ioffst(adr):(ioffst(adr) + icnt(adr) - 1));
    if(adr>1)
        sourcesNearby_left = [isradr(ioffst(adr-1):(ioffst(adr-1) + icnt(adr-1) - 1))];
    end
    if(adr<nboxes)
        sourcesNearby_right = [isradr(ioffst(adr+1):(ioffst(adr+1) + icnt(adr+1) - 1))];
    end
    sourcesNearby = [sourcesNearby_center; sourcesNearby_left; sourcesNearby_right];
    for l = sourcesNearby'
        % prepare the indices of the sparse matrix
        lInd = ((l-1)*timeLevels + 1):((l-1)*timeLevels + timeLevels);
        Idx = ((cnt - 1)*timeLevels+1):((cnt - 1)*timeLevels+timeLevels);

        if(l~=j)
            % local evaluation from other sources
            L = abs(s(j) - s(l));
            if(L<delta && L>dt)
                % get GL nodes and weights for the integral involving other
                % neighboring sources, but does not include the density at
                % the new step
                [otherLocal2_glnodes,otherLocal2_glwts] = glwt(0,(delta - L),gl);
                otherLocal2_interpWeights_temp = zeros(P,order);
                otherInterpNodesIdx_temp = zeros(P,order);
                for kk = 1:P % for each GL node calculate interpolated value at GL node
                    [otherInterpNodes,otherInterpNodesIdx_temp(kk,:)] = get_InterpNodes(otherLocal2_glnodes(kk), dt, order,(delta - L));
                    otherLocal2_interpWeights_temp(kk,:) = barycentricInterp(otherInterpNodes, order, otherLocal2_glnodes(kk));
                end
                % prepare the values of phi to evaluate the local
                % contribution from each of the other sources
                [eta,~] = glwt(L,delta,gl);
                phiL_to_delta = generalwindow((eta./delta),winData);

                % multiply the interpolation weights with the blending
                % function
                otherLocal2_interpWeights_temp = otherLocal2_interpWeights_temp.*(1 - phiL_to_delta(end:-1:1)).*otherLocal2_glwts;

                % shift the indices;
                otherInterpNodesIdx_temp = otherInterpNodesIdx_temp + snIdxNow + 1 - W;

                [otherLocal2_interpWeights] ...
                    = mergeWeights(otherInterpNodesIdx_temp,otherLocal2_interpWeights_temp,timeLevels);

                i_sp(Idx) = j;
                j_sp(Idx) = lInd;
                a_sp(Idx) = otherLocal2_interpWeights;

                cnt = cnt + 1 ;
                 deltaInteractions(j) = deltaInteractions(j) + 1; 

            elseif(L<dt)

                otherLocal1_wts = get_lmmwts(tn(1:order), order, order,(dt - L),dt,delta,winData,0);
                coef_other = (-beta(j)/2)*otherLocal1_wts(end); % coefficient of the sigma at (n+1)
                implicitMat(j,l) = coef_other;

                i_sp(Idx) = j;
                j_sp(Idx) = lInd;
                a_sp(Idx) = selfLocal2_interpWeights;

                % Adding the local 1 contribution
                endTerm = Idx(end)-1; % Index refering tothe current step
                a_sp((endTerm+1-(order-1)):endTerm) = a_sp((endTerm+1-(order-1)):endTerm) + otherLocal1_wts((end-(order-1)):end-1);
                cnt = cnt + 1 ;

                   dtInteractions(j) = dtInteractions(j) + 1;
            end
        elseif(j == l)
            % Adding the self local 2 contributions
            i_sp(Idx) = j;
            j_sp(Idx) = lInd;
            a_sp(Idx) = selfLocal2_interpWeights;

            % Adding the self local 1 contributions
            endTerm = Idx(end)-1; % Index refering to the current step
            a_sp((endTerm+1-(order-1)):endTerm) = a_sp((endTerm+1-(order-1)):endTerm) + selfLocal1_wts((end-(order-1)):end-1);
            cnt = cnt + 1 ;
        end

    end % end of l loop
end % end of j loop

% decompose the matrix to speed up solution in the time-stepping loop
implicitMat = decomposition(implicitMat);

A_sp = sparse(i_sp,j_sp,a_sp,M,(M*timeLevels));

clear i_sp j_sp a_sp;

end % end of function