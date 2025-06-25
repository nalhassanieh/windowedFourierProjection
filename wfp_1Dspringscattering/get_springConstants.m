function beta = get_springConstants(betaMax,M,uniform_beta)

if(uniform_beta == 1)
    beta = betaMax*ones(M,1);
else
    beta = 0.1 + (betaMax-0.1)*rand(M,1); % constant for the spring
end

end