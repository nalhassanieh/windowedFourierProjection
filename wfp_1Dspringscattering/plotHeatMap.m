function plotHeatMap(u,ue,x,tn,s,figFile,order,tFinal,prechar,savePlot,solnType,logScale,dataFile)
% PLOTHEATMAP plots a heatmap of a given solution u(x,t)
% 
% plotHeatMap(u,ue,x,tn,figFile,order,tFinal,prechar,savePlot)
%  returns a heatmap of u with respect to x and t given discrete 'x' grid
%  and 'tn' time grid. 'figFile' is the file where the figure will be
%  saved, 'order' and 'tFinal' are needed to record the order of accuracy 
%  and the final time in the file name. 'prechar' is a character used to
%  precede the file name if needed for special cases. 'savePlot' = 1 to
%  save plot in the given figFile. 'solnType' type of testing solution used
%  ms (manufactured solution), true (true scattering example).

M = length(s);

colormap(jet(256));

figure(1);
if(strcmp(solnType,'ms'))
    if(logScale == 1)
        U = log10(abs(u))'; 
        climInt = [-12 0];
        solTitle = 'log_{10}|u(x,t)|';
    else 
        U = abs(u)';
        climInt = [0 0.5];
        solTitle = '|u(x,t)|';
    end
    imagesc(x,tn,U);  hold on;
    plot(s,0,'|r','MarkerSize',10); hold off;
    axis xy; clim(climInt);
    xlabel x; ylabel t; colorbar; title(sprintf('M = %d, order = %d, $%s$',M,order,solTitle));
else % in this case ue = u_incident 
    if(logScale == 1)
        U = log10(abs(u + ue))';
        climInt = [-12 0]; 
        solTitle = 'log_{10}|u_{tot}(x,t)|';
    else 
        % U = abs(u + ue)';
        % climInt = [0 0.05];
        % solTitle = '|u_{tot}(x,t)|';

                U = (u + ue)';
        climInt = [-0.005 0.005];
        solTitle = 'u_{tot}(x,t)';
    end
    imagesc(x,tn,U);  hold on;
    plot(s,0,'|r','MarkerSize',10);
    % xline(s,'--r');
    hold off;
    axis xy; clim(climInt);
    xlabel x; ylabel t; colorbar; title(sprintf('M = %d, order = %d, $%s$',M,order,solTitle));
end

if(savePlot == 1)
    fileName = sprintf('%s/%sDSheatmap%d_%d_%d',figFile,prechar,order,round(tFinal),M);
    print(fileName,'-depsc');

    fileName2 = sprintf('%s/%sDSheatmapDATA%d_%d_%d',dataFile,prechar,order,round(tFinal),M);
    save(fileName2);
end

% plot the heatmap of the error
if(strcmp(solnType,'ms'))
    figure(2);
    imagesc(x,tn,log10(abs(u-ue))'); axis xy; clim([-12 0]); hold on;
    plot(s,0,'|r','MarkerSize',10); hold off;
    xlabel x; ylabel t; colorbar; title(sprintf('M = %d, order = %d, log_{10}|u(x,t)-u_e(x,t)|',M,order));

    if(savePlot == 1)
        fileName = sprintf('%s/%sDSheatmap%d_%d_%d',figFile,prechar,order,round(tFinal),M);
        print(fileName,'-depsc');
    end
end

end