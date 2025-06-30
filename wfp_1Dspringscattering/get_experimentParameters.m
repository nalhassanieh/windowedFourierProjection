function [plotOpt,saveWorkspaceOpt,addErrResults,savePlot,...
    logScale,printTimeStep,evalSol,order,numResolutions,solnType,tFinal,...
    Nx,Nt_sol,M,maxNumNeighbors,src_dmn,ds,...
    uniform_sgrid,betaMax,mu,dt0_min,uniform_beta,keepdt]...
    =get_experimentParameters(expNum)

switch expNum
    case 1

        % testing buttons
        plotOpt = 0;        % '1' to plot sol at final time, '2' heat maps
        saveWorkspaceOpt = 0;
        addErrResults = 0;
        savePlot = 0;
        logScale = 0;          % '1' to use log scale in heat maps
        printTimeStep = 0;        % '1' to print out each time step
        evalSol = 1;

        % parameters
        order = 2;            % order of accuracy = interpolation
        numResolutions = 3;   % num of grid resolutions
        solnType = 'ms';      % type of solution 'ms' or 'true'

        tFinal = 3*pi;        % final time

        %%% Choose number of spatial and temporal set for solution computation
        Nx = 10;              % size of spatial grid for testing and plotting
        Nt_sol = 10;          % size of the time grid for testing and plotting

        %%% sources
        M = 2;               % number of sources
        maxNumNeighbors = NaN; % set max number of neighbors
        src_dmn = [-1,1];
        ds = 1e-4;          % set min distance between sources
        uniform_sgrid = 1;

        %%% spring constants
        betaMax = 3; uniform_beta = 1;

        %%% incident pulse parameter
        mu = 50;
        keepdt = 1;

        dt0_min = 0.1;
    case 2

        % testing buttons
        plotOpt = 0;        % '1' to plot sol at final time, '2' heat maps
        saveWorkspaceOpt = 1;
        addErrResults = 0;
        savePlot = 0;
        logScale = 0;          % '1' to use log scale in heat maps
        printTimeStep = 1;        % '1' to print out each time step
        evalSol = 1;

        % parameters
        order = 8;            % order of accuracy = interpolation
        numResolutions = 2;   % num of grid resolutions
        solnType = 'true';      % type of solution 'ms' or 'true'

        tFinal = 10*pi;        % final time

        %%% Choose number of spatial and temporal set for solution computation
        Nx = 400;              % size of spatial grid for testing and plotting
        Nt_sol = 400;          % size of the time grid for testing and plotting

        %%% sources
        M = 150;               % number of sources
        maxNumNeighbors = 10; % set max number of neighbors
        src_dmn = [-2,2];
        ds = 1e-4;          % set min distance between sources
        uniform_sgrid = 0;

        %%% spring constants
        betaMax = 10; uniform_beta = 0;

        %%% incident pulse parameter
        mu = 30;
        keepdt = 0;

        dt0_min = 0.1;

    case 3

        % testing buttons
        plotOpt = 0;        % '1' to plot sol at final time, '2' heat maps
        saveWorkspaceOpt = 1;
        addErrResults = 0;
        savePlot = 0;
        logScale = 0;          % '1' to use log scale in heat maps
        printTimeStep = 1;        % '1' to print out each time step
        evalSol = 1;

        % parameters
        order = 8;            % order of accuracy = interpolation
        numResolutions = 2;   % num of grid resolutions
        solnType = 'true';      % type of solution 'ms' or 'true'

        tFinal = 3*pi;        % final time

        %%% Choose number of spatial and temporal set for solution computation
        Nx = 800;              % size of spatial grid for testing and plotting
        Nt_sol = 800;          % size of the time grid for testing and plotting

        %%% sources
        M = 1000;               % number of sources
        maxNumNeighbors = 50; % set max number of neighbors
        src_dmn = [-2,2];
        ds = 1e-4;          % set min distance between sources
        uniform_sgrid = 0;

        %%% spring constants
        betaMax = 3; uniform_beta = 0;

        %%% incident pulse parameter
        mu = M^2;
        keepdt = 0;

        dt0_min = 0.1;
    case 4
        % testing buttons
        plotOpt = 0;        % '1' to plot sol at final time, '2' heat maps
        saveWorkspaceOpt = 1;
        addErrResults = 0;
        savePlot = 0;
        logScale = 0;          % '1' to use log scale in heat maps
        printTimeStep = 0;        % '1' to print out each time step
        evalSol = 1;

        % parameters
        order = 8;            % order of accuracy = interpolation
        numResolutions = 2;   % num of grid resolutions
        solnType = 'true';      % type of solution 'ms' or 'true'

        tFinal = 3*pi;        % final time

        %%% Choose number of spatial and temporal set for solution computation
        Nx = 800;              % size of spatial grid for testing and plotting
        Nt_sol = 800;          % size of the time grid for testing and plotting

        %%% sources
        M = 10000;               % number of sources
        maxNumNeighbors = 50; % set max number of neighbors
        src_dmn = [-2,2];
        ds = 1e-4;          % set min distance between sources
        uniform_sgrid = 0;

        %%% spring constants
        betaMax = 3; uniform_beta = 0;

        %%% incident pulse parameter
        mu = M^2;
        keepdt = 0;

        dt0_min = 0.1;
    case 5 % FP
        % testing buttons
        plotOpt = 0;        % '1' to plot sol at final time, '2' heat maps
        saveWorkspaceOpt = 1;
        addErrResults = 0;
        savePlot = 0;
        logScale = 0;          % '1' to use log scale in heat maps
        printTimeStep = 0;        % '1' to print out each time step
        evalSol = 1;

        % parameters
        order = 8;            % order of accuracy = interpolation
        numResolutions = 2;   % num of grid resolutions
        solnType = 'true';      % type of solution 'ms' or 'true'

        tFinal = 30*pi;        % final time

        %%% Choose number of spatial and temporal set for solution computation
        Nx = 400;              % size of spatial grid for testing and plotting
        Nt_sol = 500;          % size of the time grid for testing and plotting

        %%% sources
        M = 2;               % number of sources
        maxNumNeighbors = 2; % set max number of neighbors
        src_dmn = [-0.5,0.5];
        ds = 1e-4;          % set min distance between sources
        uniform_sgrid = 1;

        %%% spring constants
        betaMax = 100; uniform_beta = 1;

        %%% incident pulse parameter
        mu = 5;
        keepdt = 0;

        dt0_min = 0.02; % RBC constraint
        dt0_min = min(0.01,dt0_min); % stability constraint

    case 6 % FP M = 10
        % testing buttons
        plotOpt = 0;        % '1' to plot sol at final time, '2' heat maps
        saveWorkspaceOpt = 0;
        addErrResults = 0;
        savePlot = 0;
        logScale = 0;          % '1' to use log scale in heat maps
        printTimeStep = 1;        % '1' to print out each time step
        evalSol = 1;

        % parameters
        order = 8;            % order of accuracy = interpolation
        numResolutions = 3;   % num of grid resolutions
        solnType = 'true';      % type of solution 'ms' or 'true'

        tFinal = 3*pi;        % final time

        %%% Choose number of spatial and temporal set for solution computation
        Nx = 10;              % size of spatial grid for testing and plotting
        Nt_sol = 10;          % size of the time grid for testing and plotting

        %%% sources
        M = 2;               % number of sources
        maxNumNeighbors = 10; % set max number of neighbors
        src_dmn = [-2,2];
        ds = 1e-4;          % set min distance between sources
        uniform_sgrid = 1;

        %%% spring constants
        betaMax = 100; uniform_beta = 1;

        %%% incident pulse parameter
        mu = 5;
        keepdt = 0;

        dt0_min = 0.02; % RBC constraint
        dt0_min = min(0.001,dt0_min); % stability constraint
    case 7 % FP M = 200
        % testing buttons
        plotOpt = 0;        % '1' to plot sol at final time, '2' heat maps
        saveWorkspaceOpt = 1;
        addErrResults = 0;
        savePlot = 0;
        logScale = 0;          % '1' to use log scale in heat maps
        printTimeStep = 0;        % '1' to print out each time step
        evalSol = 1;

        % parameters
        order = 8;            % order of accuracy = interpolation
        numResolutions = 2;   % num of grid resolutions
        solnType = 'true';      % type of solution 'ms' or 'true'

        tFinal = 40*pi;        % final time

        %%% Choose number of spatial and temporal set for solution computation
        Nx = 800;              % size of spatial grid for testing and plotting
        Nt_sol = 800;          % size of the time grid for testing and plotting

        %%% sources
        M = 200;               % number of sources
        maxNumNeighbors = 20; % set max number of neighbors
        src_dmn = [-2,2];
        ds = 1e-4;          % set min distance between sources
        uniform_sgrid = 1;

        %%% spring constants
        betaMax = 1; uniform_beta = 1;

        %%% incident pulse parameter
        mu = 5;
        keepdt = 0;

        dt0_min = 0.02; % RBC constraint
        dt0_min = min(0.01,dt0_min); % stability constraint
    otherwise
        disp('experiment number inapplicable');
end

end

