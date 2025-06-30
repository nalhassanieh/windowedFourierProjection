%% determine the file name 
if (uniform_sgrid == 1)
    str = 'uni'; else; str = 'rand'; end

% Save only the selected variables
fileName = sprintf('%s/WFPworkspaceO%d_T%d_%sM%d_b%d',dataFile,order,round(tFinal),str,M,betaMax);
    matlabDataFile =strcat(fileName,'.mat');
    
if (saveWorkspaceOpt == 1)
    % all variables in the workspace
    vars = who;

    % variables to exclude
    exclude = {'A_sp','It','RHS','deltaInteractions',...
        'dtInteractions','evalSol','gl','gn','historySum',...
        'implicitMat','interpAndGLWeights_sol','interpNodesIdx_sol',...
        'interpShift_sol','it1','it2','p','p0','phi_tn','phit_tn',...
        'printTimeStep','q','q0','saveWorkspaceOpt',...
        'sn_v','timeSteppingTime','uniform_beta','uniform_sgrid',...
        'winData'};

    % Find variables to save
    vars_to_save = setdiff(vars, exclude);
    save(fileName, vars_to_save{:});
    fprintf('saved workspace to %s\n',matlabDataFile);

elseif(addErrResults == 1)
    if(isfile(matlabDataFile))
        save(fileName, 'sc_err', '-append');
        fprintf('added error to the data file in %s\n',matlabDataFile);
    else
        fprintf('Did not add the error to the data file %s\n',matlabDataFile);
    end
end


