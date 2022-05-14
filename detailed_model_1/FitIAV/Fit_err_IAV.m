function err = Fit_err_IAV(data_h, data_m, data_neu, data_nk, data_v, data_t, data_te, data_il6, data_ccl2, parFit, parBase)
    
    % meaning of inputs
    % data*(1): time; data*(2): value of ZT23; data*(3): value of ZT11
    % parFit: 74 parameters
    
    % translate parFit to par in the odes
    log_par_ind = [1:40 49:59];
    par_IAV = [parFit; parBase];
    for i = log_par_ind
        par_IAV(i) = 10 .^ parFit(i);
    end
         
    % solve odes
    tmax = 250;
    tspan = 0:1:tmax;
    y0 = zeros(13, 1);
    y0(4) = 40;
    y0(1) = parBase(1);
    y0(5) = parBase(2);
    y0(6) = parBase(3);
    y0(11) = parBase(4);
    y0(12) = 10;
    
    
    % For ZT23
    [t, y_23] = ode15s(@ODE_IAV, tspan, y0, [], par_IAV);

    y_23 = real(y_23);

    err = 0;
    
    % cost function of V - log error
    for i = 1:8
        % err of ZT23
        ind_23 = t == data_v(1, i);
        err = err + ((Safe_log10(y_23(ind_23, 4)) - Safe_log10(data_v(2, i))) ...
            / max(Safe_log10(data_v(2, :)))) .^ 2;
    end

    % cost function of H, M, Neu and NK - non log error
    for i = 1:4
        % err of ZT23
        ind_23 = t == data_m(1, i);
        err = err + ((y_23(ind_23, 1) - data_h(2, i)) ...
            / max(data_h(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 5) - data_m(2, i)) ...
            / max(data_m(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 6) - data_neu(2, i)) ...
            / max(data_neu(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 11) - data_nk(2, i)) ...
            / max(data_nk(2, :))) .^ 2;
    end
    
    % cost function of T and TE - non log error
    for i = 1:6
        % err of ZT23
        ind_23 = t == data_t(1, i);
%         err = err + ((y_23(ind_23, 12) - data_t(2, i)) ...
%             / max(data_t(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 13) - data_te(2, i)) ...
            / max(data_te(2, :))) .^ 2;
    end
    
    % cost function of IL6 and CCL2 - non log error
    for i = 1:4
        % err of ZT23
        ind_23 = t == data_il6(1, i);
        err = err + ((y_23(ind_23, 7) - data_il6(2, i)) ...
            / max(data_il6(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 9) - data_ccl2(2, i)) ...
            / max(data_ccl2(2, :))) .^ 2;
    end

end