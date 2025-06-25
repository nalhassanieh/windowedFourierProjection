startMatlabFile

%% Set file directories to save data, figures, tables and movies
addpath ./utils;
dataFile = './data';
figFile = './fig';

%% Testing Buttons
plotOpt = 0;        % '1' to plot sol at final time, '2' heat maps
saveWorkspaceOpt = 0; 
addErrResults = 0; 
savePlot = 0; 
logScale = 0;          % '1' to use log scale in heat maps   
printTimeStep = 0;        % '1' to print out each time step
evalSol = 1; 

%% Problem parameters
%%% order and test solution
order = 2;            % order of accuracy = interpolation 
numResolutions = 3;   % num of grid resolutions
solnType = 'ms';      % type of solution 'ms' or 'true'

%%% domain and final time
ax = -pi; bx = pi;    % Domain
tFinal = 3*pi;        % final time
tsStop = -1;          % stop time stepping at tsStop. ('-1' continue to tFinal)

%%% Choose number of spatial and temporal set for solution computation
Nx = 10;              % size of spatial grid for testing and plotting
Nt_sol = 10;          % size of the time grid for testing and plotting
tSOL = tFinal/Nt_sol; % time-step for the solution evaluation

%%% sources
M = 10;               % number of sources
maxNumNeighbors = 10; % set max number of neighbors
src_dmn = [-1,1]; 
ds = 1e-4;          % set min distance between sources
uniform_sgrid = 0; 
s = get_sourceGrid(src_dmn,ds,M,uniform_sgrid); 

%%% spring constants
betaMax = 3; uniform_beta = 0; 
beta = get_springConstants(betaMax,M,uniform_beta); 

%%% incident pulse parameter
mu = 30;

%%% window
tol = 1e-12;    % error tolerance
[winData,W] = setup_generalwindow(tol); 

%% Compute some values for the window function and GL
P = W; % GL nodes 
gl = glwt_prep(P); % GL nodes and weights on [-1,1]

%% Manufactured density values for testing or incident wave for true scattering 
if(strcmp(solnType,'ms'));[h0,dataParam] = manufacturedSolution(M,tFinal);
else; [h0,dataParam] = get_uIncidentInfo(mu); end 

% fix initial time step 
[h0,typNumOfNeighbors] = fixDtBasedOnTypNeighborsNum(maxNumNeighbors,src_dmn,W,M,h0);
Nt0 = ceil(tFinal/h0); if(mod(Nt0,2) == 1); Nt0 = Nt0 + 1; end

% get the spatial grid
x = get_xgrid(ax,bx,src_dmn,Nx,h0,W);

%% Print title
printTitle(tFinal, tol, W, P, M, typNumOfNeighbors, order, solnType);
prechar = sprintf('%s%d_%d_mn%d_',solnType,betaMax,dataParam.mu(1),maxNumNeighbors); % char/string to precede title of fig or tab

%% Grid resolution study
if(evalSol == 1)
    err = zeros(1,numResolutions);
    h   = zeros(1,numResolutions);
    ug  = cell(1,numResolutions);
    t_cell = cell(1,numResolutions); 
    x_cell = cell(1,numResolutions); 
    sn_cell = cell(1,numResolutions); 
end

for m = 1:numResolutions
    N = Nt0*(2^(m-1));

    dt = 2*pi/N;
    Nt = ceil(tFinal/dt);
    dt = tFinal/Nt;          % timestep
    h_RBC = 2*pi/N;

    tgh = W + order;
    tInitial = -tgh*dt;      % initial time
    Ntg = Nt + 1 + tgh;      % total number of time-steps
    tn = linspace(tInitial ,tFinal, Ntg)'; % time grid
    it1 = tgh + 1; it2 = Ntg; % interior time indices (start from zero)
    It = it1:it2;            % time indices starting from zero
    delta = dt*W;

    %%% time parameter for self convergence study
    if(m == 1); fixSolCnt = ceil((Nt + 1)./Nt_sol);
        else; fixSolCnt = 2*fixSolCnt; end
        Nt_sol = ceil((Ntg-tgh)/fixSolCnt);

    % frequency paramenters
    K0 = N/2;                % max wave number
    K = [0:(K0-1),(-K0:-1)]; % frequency set
    zeroLoc = 1;             % where k = 0 is located

    % set up tRBC
    tRBC = W*h_RBC; % choose tFinal to be a multiple of pi

    %% Allocate space for time-stepping
    % let sn hold the density function values at any time tn
    timeLevels = tgh + 2; 
    sn = zeros(timeLevels,M);% store tgh + current time + next time
    snIdxNow = tgh+1;     % index of sn at the current time

    % Let an hold the values of the Fourier coefficients, and bn
    % correspond to time derivative of an, needed in the history treatment
    an = zeros(N,1);
    bn = an; anp1 = an; bnp1 = bn;

    % Allocate space for computed solution for heat maps
    u = zeros(Nx,Nt_sol); % u(x,t) to generate heat map
    ue = zeros(Nx,Nt_sol);
    tn_sol = zeros(1,Nt_sol);
    sn_sol = zeros(M,Nt_sol); 

    %% Prepare for time-stepping
    t = dt; % this represents the current time

    tic
    % prepare nodes and weights needed for the evaluation of self local
    % integrals and other local integrals
    [implicitMat,A_sp,deltaInteractions,dtInteractions] = ...
        prepNodesAndWeights(tn(It),dt,beta,gl,s,src_dmn,snIdxNow,P,M,W,order,...
         timeLevels,typNumOfNeighbors,winData);

    % prepare the nodes and weights for the solution local evaluation
    if(evalSol == 1)
        [interpNodesIdx_sol,interpAndGLWeights_sol,interpShift_sol] = prepForLocalSolEval(x,dt,gl,s,M,P,W,order,winData);
    end

    % prepare integrals needed in the evaluation of Fourier coefficients
    [p,q,p0,q0,sn_hat,phi_tn,phit_tn] = prepValuesForFourierCoefs(sn,snIdxNow,dt,K,N,W,gl,s,tol,winData);

    tm.prep = toc;
    %% Time Stepping
    tRBC_now = 0; tSOL_now = tSOL; n_sol = 1; totalSteps = 1;
    tic
    for n = It(2:(end-1))  % note: t(n) = tInitial + (n-1)*dt = (-tgh + (n-1))*dt

        % get the history contribution
        [anp1,bnp1,historySum,tRBC_now] = get_historyContribution(t,an,bn,sn,zeroLoc,K,...
            dt, tRBC, W, snIdxNow, p,q,p0,q0,sn_hat,tRBC_now,phi_tn,phit_tn,s,tol);

        % Get the gn for the spring scattering conditions
        gn = get_g(t,dt,s,dataParam,M,beta,solnType);

        % Perform a sparse matvec to evaluate integrals
        sn_v = reshape(sn,[],1);
        RHS = gn + (beta/2).*((A_sp*sn_v)+ (1/pi)*historySum);

        % solve the implicit system
        sn(end,:) = implicitMat\RHS;

        % Take the real part of the density
        sn = real(sn);

        % new time
        t = tn(n+1);

        if(evalSol == 1)
            [u,ue,tn_sol,n_sol,tSOL_now,sn_sol] = get_solnAndErr(u,ue,tn_sol,n_sol,...
                tSOL_now,x,t,snIdxNow,sn,anp1,interpNodesIdx_sol,...
                interpAndGLWeights_sol,interpShift_sol,s,M,delta,tol,...
                tSOL,solnType,dataParam,fixSolCnt,n,tgh,sn_sol);
        end

        % update values for the next time step
        [an,bn,sn,sn_hat] = prepForNextStep(anp1,bnp1,sn,sn_hat,snIdxNow,tgh,s,N,tol);

        if(printTimeStep == 1)
            fprintf('t = %1.4e\n',t);
        end

        if(totalSteps == tsStop)
            break;
        end
        totalSteps = totalSteps + 1;
    end

    timeSteppingTime = toc;
    tm.ts = timeSteppingTime/totalSteps;

    %% Error Analysis
    if(evalSol == 1)
        u = real(u);
        ug{m} = u; 
        h(m) = dt;
        t_cell{m} = tn_sol; 
        x_cell{m} = x'; 
        sn_cell{m} = sn_sol; 

        % Check the error at all time
        if(strcmp(solnType,'ms'))
            err(m) = max(max(abs(u-ue))); % max-norm error
            rate = printOutResults(m,t,Nt,K0,dt,err,h);
        end

        % Plot and save into a file if needed
        if(plotOpt == 1 && strcmp(solnType,'ms'))
            plotSolnAndError(t,x,u(:,end),ue,figFile,prechar,order,tFinal,savePlot);
            pause(.1);
        elseif(plotOpt == 2)
            plotHeatMap(u,ue,x,tn_sol,s,figFile,order,tFinal,prechar,savePlot,solnType,logScale,dataFile);
            pause(1);
        end
    end
end

checkMatchingGrids(t_cell,x_cell,numResolutions);

if(numResolutions>2)
    [sc_err,sc_rate] = selfConvergence(ug,h,numResolutions);
elseif(numResolutions==2)
    sc_err = maxNorm(ug{2} - ug{1});
    fprintf('self convergence error = %1.2e\n',sc_err);
end

if(saveWorkspaceOpt == 1 || addErrResults == 1)
    saveWorkspace;
end