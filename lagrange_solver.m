function [ u v ] = lagrange_solver( )
%lagrange solver: solves a linear partial differential eqution of the 
%form Pp + Qq + R = 0 by lagrange's method

syms x y z p q dx dy c1 c2
%the equation to be solved
lhs = (y^2 *z /x) * p + z*x*q;
rhs = y^2;

C = coeffs(lhs, [q p]);
P = C(1);
Q = C(2);
R = rhs;

%seperate variables
P = P*x / z;
Q = Q*x / z;

%Integrating
u = int(P, y) == int(Q,x) + c1;
%consider first and last fractions
P = C(1);
P = y^2 * z / P;
R = y^2 * z / R;
v = int(P,x) == int(R,z) + c2;
end