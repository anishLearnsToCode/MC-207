function [sols, vals] = linear_IVP_system_solver(A, B, C)
    syms t lambda
    n = length(A);
    [V, D] = eig(A);
    eigenValues = diag(D);
    constants = reshape(sym('c%d', 1:n), n, 1);
    unique_EigenValues = unique(eigenValues);
    mults = histc(eigenValues, unique_EigenValues);
    sols = sym('x%d', [1 n]);
    
    if(length(unique_EigenValues) ~= length(eigenValues))
        %For representing Eigen Values
        i = 1;
        ch_Mat = A - lambda * eye(n);
        V = vpa(V);
        
        while(i <= n)
            [pos] = find(unique_EigenValues == eigenValues(i));
            if(mults(pos) > 1)
                EigenVector = V( : , i);
                a_Mat = subs(ch_Mat, eigenValues(i));
                
                for j = 1 : mults(pos)
                    V(:, i) = V(:, i) .* (t ^ (j-1));
                    P = inv(a_Mat ^ (j-1)) * EigenVector;
                    V(:, i) = V(:, i) + P;
                    i = i + 1;
                end
            else
                i = 1+1;
            end
        end
    end
    
    for i = 1 : n
        sols(i) = (V(i, :) .* exp(eigenValues' * t)) * constants;
    end
    
    vals = solve(subs(sols, t, B) == C);
    constantNames = fieldnames(vals);
    
    %The final solution 
    for i = i : n
        sols = subs(sols, constants(i), vals.(constantNames{i}));
    end