function gn = get_g(t,dt,s,dataParam,M,beta,solnType)

if(strcmp(solnType,'ms'))
    U = get_ExactSol(s,(t + dt),dataParam.mu,dataParam.t0,M,s);  % get the manufactured solution at xval
    gn = -dataParam.sig(t + dt) - beta.*U;
else
    uin = get_uIncident(s,(t + dt),dataParam.mu,dataParam.t0);
    gn = beta.*uin;
end

end