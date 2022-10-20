function y = solve_odes(par_base, par_consider_idx, par_consider, y0, tspan)
    
    % replace parameters considered with par_consider
    tmp = 1;
    for i = par_consider_idx
        par_base(i) = par_consider(tmp);
        tmp = tmp + 1;
    end
    
    % solve odes
    [t, y] = ode15s(@ODE_IAV, tspan, y0, [], par_base);
    y = real(y);
   
end