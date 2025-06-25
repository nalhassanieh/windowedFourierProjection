function checkMatchingGrids(t_cell,x_cell,numResolutions)

tGridErr = 0; xGridErr = 0; 
for m = 2:numResolutions
    tGridErr = max(tGridErr,maxNorm(t_cell{m} - t_cell{1}));
    xGridErr = max(xGridErr,maxNorm(x_cell{m} - x_cell{1}));
end

tol = 1e-12;
if(tGridErr>tol || xGridErr>tol)
    warning('Time or space grids do not match to report self-convergence');
end

end
