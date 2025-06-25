function [interpNodes,interpNodesIdx] = get_InterpNodes(auxNode, dt, order, endTime)
% GETINTERPNODES retrieves the interpolation nodes on the time-grid to
%  interpolate unknown function values at auxiliary nodes
%
% [interpNodes,interpNodesIdx] = get_InterpNodes(auxNode, dt, order, endTime)
%   returns the interpolation nodes on the time-grid and their
%   corresponding indices to interpolate function values at a given
%   auxiliary node (off the time-grid). The function ensures that
%   interpolation nodes do not go beyong a given endTime
% 
% Inputs: 
%  auxNode: given value that may not coincide with the uniform time-grid
%  dt: time-step
%  order: equals the number of interpolation nodes used
%  endTime: ensure interpolation nodes < endTime
%
% Outputs: 
%  interpNodes: a vector of size order that carries the interpolation nodes
%     needed to determine function values at auxNode
%  interpNodesIdx: a vector of size order that carries the interpolation
%     nodes indices on the time-grid


%%% This is the one-sided treatment (towards the left)
% rightGridPoint_idx = ceil(glnode/dt);
% rightGridPoint = rightGridPoint_idx*dt; % find the nearest node on the uniform grid
% interpNodes = linspace(rightGridPoint  - (order - 1)*dt,rightGridPoint,order)';
% interpNodesIdx = linspace((rightGridPoint_idx - (order - 1)),rightGridPoint_idx,order)';

if nargin==0, test_getInterpNodes; return; end

%%% This is the centered/skewed treatment
m = order; % set m to equal the order

tau_idx = round(auxNode/dt); % index of nearest node on time-grid
tau = tau_idx*dt;            % nearest node on time-grid
endTime_idx = round(endTime/dt); % index of endTime
k = endTime_idx - tau_idx;       % how far is nearest node from endTime

if(mod(m,2) == 0) % if m is even
    if(k>=((m-2)/2)) % centered interpolation treatment 

        if(auxNode>tau)
            % rightward skew: add interpolation point to the right end
            startingNode = tau - ((m-2)/2)*dt;
            startingIdx = tau_idx - ((m-2)/2);

            % if last node is outside domain: add interpolation point to
            % the left end
            if((tau+((m/2)*dt))>endTime) 
                startingNode = tau - (m/2)*dt;
                startingIdx = tau_idx - (m/2);
            end

        else % leftward skew: add interpolation point to the left end
            startingNode = tau - (m/2)*dt;
            startingIdx  = tau_idx - (m/2);
        end
    else 
        % skewed treatment: when centered interpolation is not possible
        % because interpolation nodes become larger than endTime
        startingNode = tau - (m-k-1)*dt;
        startingIdx  = tau_idx - (m-k-1);
    end
else % if m is odd
    if(k>=((m-1)/2)) % centered interpolation
        startingNode = (tau - ((m-1)/2)*dt);
        startingIdx = tau_idx - ((m-1)/2);
    else % skewed treatment
        startingNode = tau - (m-k-1)*dt;
        startingIdx  = tau_idx - (m-k-1);
    end
end

% interpolation Nodes: startingNode, startingNode + dt,...
%                            ...,startingNode + (m-1)*dt

interpNodes = (startingNode:dt:(startingNode + ((m-1)*dt)))';
interpNodesIdx = (startingIdx:(startingIdx + (m-1)))';

end

function test_getInterpNodes
clf;

startTime = 0; endTime = 1; 
t1 = .25; t2 = 1; 
N = 16; 
order = 5; 
saveToFile = 0; % set to '1' 
saveFile = '/Users/nalhassanieh/Desktop';

tn = linspace(startTime,endTime,N); 
dt = tn(2) - tn(1); 
glnodes = lgwt(N,t1,t2); 

% plot integration domain and GL nodes
figure(1);
hAxes = axes('NextPlot','add',...
    'DataAspectRatio',[1 1 1],...
    'XLim',[startTime endTime],...
    'YLim',[0 eps],...
    'Color','none');
p1 = plot(glnodes,0,'bo','MarkerSize',10);
p2 = plot(tn,0,'r|','MarkerSize',10);
p3 = plot([t1;t2],0,'kx','MarkerSize',10);
l_temp = [p1(1);p2(1);p3(1)];
% legend(l_temp,'GL nodes','time-grid','endpoints','interpreter','latex','location','north');
% title('GL nodes for [t_1,t_2] on time-grid'); 
set(gca,'XTick',[t1 t2], 'YTick', []);
xticklabels({'$t_1$','$t_2$'})
set(gca,'TickLabelInterpreter','latex');

if(saveToFile == 1)
    saveas(gcf,sprintf('%s/GLnodesOnTimeGrid',saveFile),'epsc');
end

pos = {'right','center','left'};

cnt = 1; 
for glnodeNum = [2,8,15]
    auxNode = glnodes(glnodeNum);
    [interpNodes,~] = get_InterpNodes(auxNode, dt, order, endTime);

    % plot one GL node and interpolation nodes
    figure;
    hAxes = axes('NextPlot','add',...
        'DataAspectRatio',[1 1 1],...
        'XLim',[startTime endTime],...
        'YLim',[0 eps],...
        'Color','none');
    p1 = plot(glnodes(glnodeNum),0,'bo','MarkerSize',10);
    p2 = plot(interpNodes,0,'*','color',[0.4660 0.6740 0.1880],'MarkerSize',10);
    p3 = plot(tn,0,'r|','MarkerSize',10);
    p4 = plot([t1;t2],0,'kx','MarkerSize',10);
    l_temp = [p1(1);p2(1);p3(1);p4(1)];
    % legend(l_temp,'GL node','interpolation nodes','time-grid','endpoints','interpreter','latex','location','north');
    % title(sprintf('%s GL node interpolation',pos{cnt})); 
    set(gca,'XTick',[t1 t2], 'YTick', []);
    xticklabels({'$t_1$','$t_2$'})
    set(gca,'TickLabelInterpreter','latex');

    if(saveToFile == 1)
        saveas(gcf,sprintf('%s/exampleGLInterp_%s',saveFile,pos{cnt}),'epsc');
    end

    pause(0.1);
    cnt = cnt + 1;
end

end