function err = Fit_err(biolu_data_ko, par_fit, par_fix)

    par_fit(3) = 10 ^ par_fit(3);
    par_scg = [par_fit(1:2); par_fix(1:3); par_fit(3); par_fix(4); par_fit(4)];
    n_E = par_scg(3);
    n_I = par_scg(4);

    % initial values
    y0 = zeros(n_E+n_I+2, 1);
    y0(1) = 1;
    
    tmax = 100;
    tspan = 0:0.1:tmax;
    
    [~, y] = ode15s(@ODE_SCG, tspan, y0, [], par_scg);
    if size(y, 1) < length(tspan)
        y = zeros(length(tspan), n_E+n_I+2);
    end

    err = 0;
    for i = 1:length(biolu_data_ko)
%         idx = tspan == round(biolu_data_ko(i, 1) * 10) / 10;
        err = err + ((y(round(biolu_data_ko(i, 1) * 10), n_E+n_I+2) - biolu_data_ko(i, 2)) ...
            / max(biolu_data_ko(:, 2))) .^ 2;
    end

end