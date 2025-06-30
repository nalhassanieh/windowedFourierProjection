function get_msOrIncidentField(M,tFinal,mu,solnType)


if(strcmp(solnType,'ms'))
    override = 0; % turn on '1' to avoid restrictions on chosen test function

mu = linspace(40,50,M)'; 
t0 = linspace(1,3,M)'; 

% choose the initial time step
dt0_vec = pi./(4*sqrt(mu*log(1/eps)));
dt0 = min(dt0_vec); 

% ensure that the forcing function is bandlimited within the computational time domain
if(any(2*t0>tFinal)&&(override == 0))
    error('Forcing inappropriate for testing. For tFinal=%1.2f\nChoose other mu or tFinal',tFinal);
end

% Define the density functions to be Gaussians
sig = @(t) exp(-mu.*((t - t0).^2));
   
dataParam.sig = sig; 
dataParam.t0 = t0; 
dataParam.mu = mu; 

else

t0 = -3; 

% choose the initial time step
dt0 = pi/(4*sqrt(mu*log(1/eps)));

dataParam.mu = mu; 
dataParam.t0 = t0; 
end 




end