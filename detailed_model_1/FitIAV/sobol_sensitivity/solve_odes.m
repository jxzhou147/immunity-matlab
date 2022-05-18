function y = solve_odes( par_IAV_ini, par_consider_ind, par_consider )
    
    % replace parameters considered with par_consider
    tmp = 1;
    for i = par_consider_ind
        par_IAV_ini(i) = par_consider(tmp);
        tmp = tmp + 1;
    end
    
    % translate parFit to par in the odes
    log_par_ind = [1:40 49:59];
    par_IAV = par_IAV_ini;
    for i = log_par_ind
        par_IAV(i) = 10 .^ par_IAV_ini(i);
    end
    
    % solve odes
    tmax = 250;
    tspan = 0:1:tmax;
    y0 = zeros(13, 1);
    y0(4) = 40;
    y0(1) = par_IAV(60);
    y0(5) = par_IAV(61);
    y0(6) = par_IAV(62);
    y0(11) = par_IAV(63);
    y0(12) = 10;
    
    [t, y] = ode15s(@ODE_IAV, tspan, y0, [], par_IAV);
    
    y = real(y);
    
end

