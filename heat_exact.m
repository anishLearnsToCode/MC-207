%%heat_exact: Compute the exact symbolic rpresentation of the solution of
%%the heat equation for a uniform bar with diffusivity k, initial
%%temperature distribution f and insulated boundaries 0 < x < a, with
%%initial temperature distribution along the rod given by f(x)

function[sol, eValue] = heat_exact(k, a)
%%using seperation of variables
syms G(t) phi(x) lambda
[eValue phi coeff] = strum_liouville(a);
phi = subs(phi, lambda, eValue);
t_ode = diff(G, t) == -k*lambda*G;
t_sol = dsolve(t_ode);
sol = t_sol * phi;
end