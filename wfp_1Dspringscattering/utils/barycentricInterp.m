function interpWeights = barycentricInterp(interpNodes, npts, xval)
% BARYCENTRICINTERP returns the barycentric interpolation weights given a
% set of interpolation nodes
% 
% interpWeights = barycentricInterp(interpNodes, npts, xval) 
%  returns barycentric interpolation weights given 'interpNodes', a vector 
%  of size 'npts' containing interpolation nodes, and 'xval', a value where 
%  the barycentric interpolating polynomial is evaluated

% if nargin == 0, test_barycentricInterp; return; end

% Precompute barycentric interpolation factors
X=repmat(interpNodes,1,npts);
w = 1./prod(-X+X.'+eye(npts),1);

tol = 1.0e-20;
xvalminusInterpNodes = xval - interpNodes';
rr = prod(xvalminusInterpNodes);

interpWeights = rr*w./(xvalminusInterpNodes);
interpWeights(abs(xval - interpNodes')<tol) = 1;

end

% function test_barycentricInterp
% clf;
% ax = 0; bx = 1;
% npts = 6; 
% interpNodes = linspace(ax,bx,npts); 
% 
% xval = (bx-ax)*rand + ax;
% interpWeights = barycentricInterp(interpNodes, npts, xval); 
% 
% bar(interpWeights); 
% title(sprintf(['Interpolation weights\n for %d equidistant interpolation...' ...
%     ' nodes\n in [%1.0f,%1.0f] evaluated at x = %1.2f'],npts,ax,bx,xval));
% 
% end