function plotSolnAndError(t,x,un,ue,figFile,prechar,order,tFinal,savePlot)
% PLOTSOLNANDERROR plots the solution and error
% 
% plotSolnAndError(t,tn,sn,src,sige,figFile,prechar,order,tFinal,plotOption)
% Inputs: 
%  tn: uniform time grid
%  sn: density grid function
%  sige: exact density function
%  figFile: file to store figures
%  prechar: character to precede name of figure

figure(1);
plot(x,un,'-o','LineWidth',2); hold on
plot(x,ue,'-x','LineWidth',2);hold off
title(sprintf('Computed vs exact at t = %1.1f',t));
axis tight;
legend('u','u_e','Location','best');
xlabel('x'); ylabel('y');
set(gca,'fontSize',16);
if(savePlot == 1)
    saveas(gcf,sprintf('%s/%sDScomputed%d_%d',figFile,prechar,order,round(tFinal)),'png');
end

figure(2);
plot(x,un-ue,'LineWidth',2);
title(sprintf('Error t = %1.1f',t));
axis tight;
legend('error','Location','best');
xlabel('t'); ylabel('y');
set(gca,'fontSize',16);
if(savePlot == 1)
    saveas(gcf,sprintf('%s/%sDSerror%d_%d',figFile,prechar,order,round(tFinal)),'png');
end

end