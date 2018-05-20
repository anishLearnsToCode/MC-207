%%Strum-Liouville: Calculate the eigenValues and eigenVectors
%%  for a strumLioville problem with biundaries [0,L] and EigenValue
%%  lambda st. X'(L)=X'(0)=0

function[eValue, eFunction, nonZero] = strum_liouville(L)
syms y(x)
syms lambda n
sprintf('solving for various conditions')
sprintf('lambda > 0')
assume(lambda > 0);

solution = dsolve(diff(y,2) + lambda*y == 0);
eFunction =solution;
diff_sol = diff(solution, x);
vals = solve(subs(diff_sol, 0)==0, subs(diff_sol, L)==0);
nonZero = vals;
sprintf('Non Zero values in the solution')
disp(vals);

%%Eigen values for this solution
eValue = [(n * pi) / L] .^ 2;
sprintf('When lambda = 0')
try
    solution = dsolve(subs(diff(y,2) + lambda*y == 0), lambda, 0);
catch
    sprintf('No non trivial solution');
end

sprintf('When lambda < 0')
assume(lambda < 0);
solution = dsolve(diff(y,2) + lambda * y == 0);
diff_sol = diff(solution, x);

%%No explicit non-trivial solutions possible
vals = solve(subs(diff_sol,0) == 0, subs(diff_sol, L) == 0)