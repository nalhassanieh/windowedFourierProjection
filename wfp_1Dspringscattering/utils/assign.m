function [ioffst, ibox, isradr,icnt] = assign(nboxes, xat, natoms)
% ASSIGN Assigns sources and targets to boxes.
%
% INPUTS:
%   nboxes - Number of boxes
%   xat - Positions of sources
%   natoms - Number of sources
%
% OUTPUTS:
%   ioffst - Offsets for sources in each box
%   ibox - Addresses of the boxes containing each target
%   center - Centers of each box
if nargin == 0, test_assign; return; end

% xat =(xat - xat(1))./(xat(end) - xat(1));

% Initialize variables
h = 1 / nboxes;
icnt = zeros(nboxes, 1);
ibox = zeros(natoms, 1);
isradr = zeros(natoms, 1);

% Find box in which each source lies and increment counter
for j = 1:natoms
    ixh = floor(xat(j)/h);
    if ixh >= nboxes
        ixh = nboxes-1;
    elseif ixh <= 0
        ixh = 0;
    end
    iadr = ixh + 1;
    icnt(iadr) = icnt(iadr) + 1;
    ibox(j) = iadr;
end

% Compute ioffst array
ioffst = zeros(nboxes, 1);
ioffst(1) = 1; 
for j = 2:nboxes
    ioffst(j) = ioffst(j-1) + icnt(j-1);
end

% Reset icnt for recording source locations in the original array
icnt(:) = 0;
for j = 1:natoms
    iadr = ibox(j);
    indx = ioffst(iadr) + icnt(iadr);
    isradr(indx) = j;
    icnt(iadr) = icnt(iadr) + 1;
end

end

function test_assign
as = -1; bs = 0; 
Ns = 10; 
s = linspace(as,bs,Ns); 
s_idx = randperm(Ns,Ns);
s  = s(s_idx);

figure(1); 
plot(s,0,'*b');

nboxes = 5; 

snew = (s - as)./(bs - as); 

[ioffst, ibox, isradr,icnt] = assign(nboxes, snew, Ns); 

plot(s,0,'b*',s(ioffst),0,'r|');

box = 3
start = ioffst(box)
isradr(ioffst(box):(ioffst(box) + icnt(box) - 1))

pause; 
end
