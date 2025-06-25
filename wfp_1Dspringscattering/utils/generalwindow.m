function [y,yp,ypp] = generalwindow(x,data)
%
% [y,yp,ypp] = generalwindow(tol,x,data)
%
% evaluates blending window on [0,1] which decays to O(tol) at 0^+,
% and 1-O(tol) at 1^-. y is vals, yp are derivs, ypp are 2nd derivs.
% yp has support in [0,1],
% x can be an array, and y, yp, ypp will be the same size as x.
%
% Note: does not return function handles!
%
% To test, call without arguments.
  
if nargin==0, test_generalwindow; return; end

tol = data.tol; 
gam = data.gam; 

% handle values x<0 and x>1...
y = 0*x; yp= 0*x; ypp = 0*x;     % also makes outputs correct sizes
y(x>=1.0) = 1.0;
% remaining x's actually need to evaluate something (not 0 or 1)...
jj = x >0.0 & x < 1.0;    % indices to eval

if isempty(data.info)    % slow evaluation
  theta = log(1/tol);    % set up window
  W = ceil(2*theta/(pi*gam));
  t = W*x(jj);           % rescaled ordinates
  [phit1,phitt1] = window(theta,W);   % get func handles for [0,W]  
  yp(jj) = W*phit1(t);           % make domain of x be [0,1]. w/ Jacobian
  ypp(jj) = W*W*phitt1(t);
  y(jj) = blending(phit1,t,tol);   % must be scalar input

else  % use data   fast eval
  t = x(jj);           % rescaled ordinates
  y(jj) = chebEval(t,data.wei,data.info);
  yp(jj) = chebEval(t,data.weip,data.info);
  ypp(jj) = chebEval(t,data.weipp,data.info);
end

end

  %%%%%%%%%%%%%
function test_generalwindow
verb = 1;
tol = 1e-12; 
% setup the data struct
data = struct('tol',tol,'gam',0.5,'info',[]); 

x = -1:1e-2:2.0;
[y,yp,ypp] = generalwindow(x,data);
if verb, 
    
    %figure; plot(x,[y;yp;ypp],'o-'); legend('y','yp','ypp'); 

    %%% testing here
    delta = 0.5; 
    t = x*delta;
    [y_new,~,~] = generalwindow((t./delta),data);
    figure; 
    plot(x,y,'o-'); xline(1); hold on; 
    plot(t,y_new,'r*'); xline(delta)
end

n=1e3; z = rand(n,1);
tic;
[Y,Yp,Ypp] = generalwindow(z,data);
fprintf("throughput of generalwindow slow = %.3g pts/sec\n",n/toc)

% test the fast cheb eval...
data = setup_generalwindow(tol);
%data
%data.info
[yf,ypf,yppf] = generalwindow(x,data);

disp("test fast vs slow eval...")
norm(y-yf,inf)
norm(yp-ypf,inf)/norm(yp,inf)
norm(ypp-yppf,inf)/norm(ypp,inf)

n=1e5; z = rand(n,1);
tic;
[Y,Yp,Ypp] = generalwindow(z,data);
fprintf("throughput of generalwindow cheb = %.3g pts/sec\n",n/toc)
end

