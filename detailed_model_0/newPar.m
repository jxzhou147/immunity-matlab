function newx = newPar(optimValues, problem)
    % generate new points using the method by Smith et al. 2018.
    
    newx = zeros(problem.nvar, 1);
    
    for i = 1:problem.nvar
        while 1
            rand_i = rand;
            g = sign(rand_i - 0.5) * optimValues.temperature(i) * ...
                ((1 + 1 / optimValues.temperature(i)) ^ abs(2 * rand_i - 1) - 1);
            newx(i) = optimValues.x(i) + g * (problem.ub(i) - problem.lb(i));
            
            if (newx(i) > problem.lb(i)) && (newx(i) < problem.ub(i))
                break;
            end
        end
    end
end