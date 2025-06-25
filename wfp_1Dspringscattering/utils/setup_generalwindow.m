function [data,W] = setup_generalwindow(tol)

% [data,W] = SETUP_GENERALWINDOW(tol)
% set up fast cheb eval for window
%
% See: generalwindow.m

gam = 0.5; data.gam = 0.5;
data.tol = tol;
data.info = {}; 

chebApproxInfo.domain = [0,1];
chebApproxInfo.ninters = max(1,ceil(0.2*log(1/tol)));  % *** hack to make sure accurate to tol!;
chebApproxInfo.nord = 16;

% table of cheb weights times function values
data.wei = mktab_wts(@(x) yfun(x,data), chebApproxInfo);
data.weip = mktab_wts(@(x) ypfun(x,data), chebApproxInfo);
data.weipp = mktab_wts(@(x) yppfun(x,data), chebApproxInfo);
data.info = chebApproxInfo;


% Some more information needed
theta = log(1/tol);    % set up window
W = ceil(2*theta/(pi*gam));

end

function y = yfun(x,data)
  [y,~,~] = generalwindow(x,data);
end
function yp = ypfun(x,data)
  [~,yp,~] = generalwindow(x,data);
end
function ypp = yppfun(x,data)
  [~,~,ypp] = generalwindow(x,data);
end
