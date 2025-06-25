function [err,rate] = selfConvergence(ug,h,M)

if nargin == 0, test_selfConvergence; return; end

if(M<3)
    error('need at least 3 resolutions to perform self-convergence study');
end

m = M-2; 

Dmp1 = ug{m+1} - ug{m};
Dmp2 = ug{m+2} - ug{m+1};

NDmp1 = maxNorm(Dmp1);
NDmp2 = maxNorm(Dmp2);
ratio = NDmp1/NDmp2;
r = h(m)/h(m+1);
rate  = log(ratio)/log(r);

C = Dmp2/abs((h(m+1)^(rate)) - (h(m+2)^(rate)));

err = zeros(1,M);
for k=1:M
    uerr = C*(h(k)^(rate));
    err(k) = maxNorm(uerr);
end

fprintf('self convergence rate = %1.2e ', rate);
fmt=['errors =' repmat(' %1.2e;',1,numel(err)) '\n'];
fprintf(fmt,err);
end

%%% test function
function test_selfConvergence

ue = @(x,t) cos(x + 3).*exp(x + 1).*log(1 + x).*sin(2*t - 2).*exp(-t);

h0 = 1e-3; 
h(1) = h0; 
h(2) = h0/2; 
h(3) = h0/4;
h(4) = h0/8; 

x = -pi:h0:pi; 
t = linspace(0,1,length(x)); 

rate = 1.5; 
numRes = 4; 
ug = cell(numRes,1); 
C = (1 + (100-1).*rand(size(x))); 
for k = 1:numRes
    ug{k} = ue(x,t) +  C*(h(k)^(rate)) +  C*(h(k)^(rate+3));
end

selfConvergence(ug,h,numRes); 

end

