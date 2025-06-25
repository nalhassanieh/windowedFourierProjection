function  s = get_sourceGrid(src_dmn,ds,M,uniform_sgrid)

as = src_dmn(1); bs = src_dmn(2); 

if(uniform_sgrid == 1)
    Ns = M;
    s_grid = linspace(as,bs,Ns);
    s = s_grid';
else
    Ns = ceil((bs - as)/ds);
    s_grid = linspace(as,bs,Ns); %s = s_grid';
    s_idx = randperm(Ns,M);
    s  = s_grid(s_idx); s = s';
end

end