%%wave: Compute the exact symbolic representation of the wave equation for
%%waves along a string of path of length a, with given initial and boundary
%%conditions

function[sol, eValue] = wave(k, a)
%%using seperaton of variables 
syms X(x) T(t) sig c1 c2 d1 d2 c3 c4

[eValue_t X coeff_X] = strum_liouville(sig);
[eValue_t T coeff_t] = strum_liouville(k*sig);
T = subs(subs(subs(T, x, t), c3, d1), c4, d2);
sol = [T * X];
end