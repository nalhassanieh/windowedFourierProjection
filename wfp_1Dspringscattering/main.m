startMatlabFile

%% Set file directories to save data, figures, tables and movies
addpath ./utils;
dataFile = './data';
figFile = './fig';

%% Problem parameters
expNum = 1; % experiment number

[plotOpt,saveWorkspaceOpt,addErrResults,savePlot,...
    logScale,printTimeStep,evalSol,order,numResolutions,solnType,tFinal,...
    Nx,Nt_sol,M,maxNumNeighbors,src_dmn,ds,...
    uniform_sgrid,betaMax,dataParam,dt0_min,uniform_beta,keepdt]...
     = get_experimentParameters(expNum); 

%%% domain 
ax = -pi; bx = pi;    % Domain
tsStop = -1;          % stop time stepping at tsStop. ('-1' continue to tFinal)

%%% sources and spring constants
s = get_sourceGrid(src_dmn,ds,M,uniform_sgrid); 
beta = get_springConstants(betaMax,M,uniform_beta); 

%%% Choose number of spatial and temporal set for solution computation
tSOL = tFinal/Nt_sol; % time-step for the solution evaluation

%%% window
tol = 1e-12;    % error tolerance
[winData,W] = setup_generalwindow(tol); 

%% Manufactured density values for testing or incident wave for true scattering 
if(strcmp(solnType,'ms'));[dt0,dataParam] = manufacturedSolution(M,tFinal);
else; [dt0,dataParam] = get_uIncidentInfo(dataParam); end 
dt0 = min(dt0_min,dt0); % domain RBC constraint 

% fix initial time step 
[dt0,typNumOfNeighbors] = fixDtBasedOnTypNeighborsNum(maxNumNeighbors,src_dmn,W,M,dt0,keepdt);
Nt0 = 2*ceil(0.5*tFinal/dt0); 
dt0 = tFinal/Nt0; 

%% Compute some values for the window function and GL*
P = W; % GL nodes 
gl = glwt_prep(P); % GL nodes and weights on [-1,1]

% get the spatial grid
x = get_xgrid(ax,bx,src_dmn,Nx,dt0,W);

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
    tm_cell = cell(1,numResolutions);
end

for m = 1:numResolutions

    dt = dt0/(2^(m-1));
    Nt = Nt0*2^(m-1);    % number of time steps

    N = 2*ceil(0.5*2*(pi/dt)); % N = number of Fourier modes case gamma = 0.5;
    % ensure N is even

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

        if(printTimeStep == 1 && mod(totalSteps,100) == 0)
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
        tm_cell{m} = tm; 

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