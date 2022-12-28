function err = Fit_err(biolu_data_ko, par_fit)

    % initial values
    y0 = zeros(4, 1);
    y0(1) = 1;
    y0(4) = 100;
    
    tmax = 100;
    tspan = 0:0.1:tmax;
    
    [~, y] = ode15s(@ODE_SCG, tspan, y0, [], par_fit);
    if size(y, 1) < length(tspan)
        y = zeros(length(tspan), 4);
    end

    err = 0;
    for i = 1:length(biolu_data_ko)
%         idx = tspan == round(biolu_data_ko(i, 1) * 10) / 10;
        err = err + ((y(round(biolu_data_ko(i, 1) * 10), 4) - biolu_data_ko(i, 2)) ...
            / max(biolu_data_ko(:, 2))) .^ 2;
    end

end