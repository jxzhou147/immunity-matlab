function err = Fit_err_clock_IAV(data_h, data_m, data_neu, data_nk, data_v, data_t, data_te, data_il6, data_ccl2, parFit, parBase)
    
    % meaning of inputs
    % data*(1): time; data*(2): value of ZT23; data*(3): value of ZT11
    % parFit: 74 parameters
    
    % translate parFit to par in the odes
    log_par_ind = [1:40 49:59 62 69 71:72];
    par_IAV = [parFit; parBase];
    for i = log_par_ind
        par_IAV(i) = 10 .^ parFit(i);
    end
    
    par_clock = load('par_clock.csv');
     
    % solve odes
    tmax = 400;
    tspan = 0:1:tmax;
    y0 = zeros(25, 1);
    y0(16) = 40;
    y0(13) = parBase(1);
    y0(17) = parBase(2);
    y0(18) = parBase(3);
    y0(23) = parBase(4);
    y0(24) = 10;
    
    
    % For ZT23
    t_IAV_23 = 115;
    [t, y_23] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_23);

    % For ZT11
    t_IAV_11 = 127;
    [t, y_11] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_11);

    y_23 = real(y_23);
    y_11 = real(y_11);

    err = 0;
    
    % cost function of V - log error
    for i = 1:8
        % err of ZT23
        ind_23 = t == (data_v(1, i) + t_IAV_23);
        err = err + ((Safe_log10(y_23(ind_23, 16)) - Safe_log10(data_v(2, i))) ...
            / max(Safe_log10(data_v(2, :)))) .^ 2;
        % err of ZT11
        ind_11 = t == (data_v(1, i) + t_IAV_11);
        err = err + ((Safe_log10(y_11(ind_11, 16)) - Safe_log10(data_v(3, i))) ...
            / max(Safe_log10(data_v(3, :)))) .^ 2;
    end

    % cost function of H, M, Neu and NK - non log error
    for i = 1:4
        % err of ZT23
        ind_23 = t == (data_m(1, i) + t_IAV_23);
        err = err + ((y_23(ind_23, 13) - data_h(2, i)) ...
            / max(data_h(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 17) - data_m(2, i)) ...
            / max(data_m(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 18) - data_neu(2, i)) ...
            / max(data_neu(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 23) - data_nk(2, i)) ...
            / max(data_nk(2, :))) .^ 2;
        % err of ZT11
        ind_11 = t == (data_m(1, i) + t_IAV_11);
        err = err + ((y_11(ind_11, 13) - data_h(3, i)) ...
            / max(data_h(3, :))) .^ 2;
        err = err + ((y_11(ind_11, 17) - data_m(3, i)) ...
            / max(data_m(3, :))) .^ 2;
        err = err + ((y_11(ind_11, 18) - data_neu(3, i)) ...
            / max(data_neu(3, :))) .^ 2;
        err = err + ((y_11(ind_11, 23) - data_nk(3, i)) ...
            / max(data_nk(3, :))) .^ 2;
    end
    
    % cost function of T and TE - non log error
    for i = 1:6
        % err of ZT23
        ind_23 = t == (data_t(1, i) + t_IAV_23);
%         err = err + ((y_23(ind_23, 24) - data_t(2, i)) ...
%             / max(data_t(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 25) - data_te(2, i)) ...
            / max(data_te(2, :))) .^ 2;
        % err of ZT11
        ind_11 = t == (data_t(1, i) + t_IAV_11);
%         err = err + ((y_11(ind_11, 24) - data_t(3, i)) ...
%             / max(data_t(3, :))) .^ 2;
        err = err + ((y_11(ind_11, 25) - data_te(3, i)) ...
            / max(data_te(3, :))) .^ 2;
    end
    
    % cost function of IL6 and CCL2 - non log error
    for i = 1:4
        % err of ZT23
        ind_23 = t == (data_il6(1, i) + t_IAV_23);
        err = err + ((y_23(ind_23, 19) - data_il6(2, i)) ...
            / max(data_il6(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 21) - data_ccl2(2, i)) ...
            / max(data_ccl2(2, :))) .^ 2;
        % err of ZT11
        ind_11 = t == (data_il6(1, i) + t_IAV_11);
        err = err + ((y_11(ind_11, 19) - data_il6(3, i)) ...
            / max(data_il6(3, :))) .^ 2;
        err = err + ((y_11(ind_11, 21) - data_ccl2(3, i)) ...
            / max(data_ccl2(3, :))) .^ 2;
    end

end