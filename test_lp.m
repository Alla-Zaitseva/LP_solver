n = 600;
m = 350;
lprog = LinearProgram(n,m);
lprog.upper_bound = +Inf;
lprog.lower_bound = -Inf;
tic
solution = solve(lprog);
fprintf("Solution found in %5e sec\n",toc)
delete(lprog)

