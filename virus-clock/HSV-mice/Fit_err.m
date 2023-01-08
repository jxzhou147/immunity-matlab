function err = Fit_err(data, par_fit)

    % initial values
    y0 = zeros(4, 1);
    y0(1) = 100000;
    y0(4) = 100;
    
    tmax = 400;
    tspan = 0:1:tmax;
    
    [~, y] = ode15s(@ODE_SCG, tspan, y0, [], par_fit);
    if size(y, 1) < length(tspan)
        y = zeros(length(tspan), 4);
    end

    err = 0;
    for i = 1:length(data)
%         idx = tspan == round(biolu_data_ko(i, 1) * 10) / 10;
        err = err + ((y(data(i, 1) + 1, 4) - data(i, 2)) ...
            / max(data(:, 2))) .^ 2;
    end

end