function [rate] = printOutResults(m,t,Nt,K0,dt,err,ht)
% PRINTOUTRESULTS prints out error results and rates in the command window
%
% printOutResults(m,t,Nt,K0,dt,err,L2err,ht)

rate = 0; 
fprintf('t=%8.2e: Nt=%4d K = %3d dt=%9.3e max-err=%8.2e',t,Nt,K0,dt,err(m));
if( m>1 )
    rate = log(err(m-1)/err(m))/log(ht(m-1)/ht(m));
    fprintf(' <strong>inf-rate=%4.2f</strong>',rate);
        fprintf('\n');
else
        fprintf('\n');
end

end