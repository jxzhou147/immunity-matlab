function err = ODE_fit_clock_IAV(data_v, data_m, data_nk, data_t, data_te, parFit)
    
    % meaning of inputs
    % data*(1): time; data*(2): value of ZT23; data*(3): value of ZT11
    % parFit: 29 parameters
    
    % translate parFit to par in the odes
    par_IAV = parFit;
    for i = [(1:12), (17:27)]
        par_IAV(i) = 10 .^ parFit(i);
    end
    
    par_clock = load('par_clock.csv');
    
    % solve odes
    tmax = 400;
    tspan = 0:1:tmax;
    y0 = zeros(19, 1);
    y0(13) = 10 ^ 7;
    y0(14) = 40;
    y0(18) = 1000;
    
    % For ZT23
    t_IAV_23 = 115;
    [t, y_23] = ode45(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_23);

    % For ZT11
    t_IAV_11 = 127;
    [t, y_11] = ode45(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV_11);

    y_23 = real(y_23);
    y_11 = real(y_11);
    
    err = 0;
    
    % cost function of V
    for i = 1:8
        % err of ZT23
        ind_23 = t == (data_v(1, i) + t_IAV_23);
        err = err + ((Safe_log10(y_23(ind_23, 15)) - Safe_log10(data_v(2, i))) ...
            / max(Safe_log10(data_v(2, :)))) .^ 2;
        % err of ZT11
        ind_11 = t == (data_v(1, i) + t_IAV_11);
        err = err + ((Safe_log10(y_11(ind_11, 15)) - Safe_log10(data_v(3, i))) ...
            / max(Safe_log10(data_v(3, :)))) .^ 2;
    end
    
%     % cost function of M and NK
%     for i = 1:4
%         % err of ZT23
%         ind_23 = t == (data_m(1, i) + t_IAV_23);
%         err = err + ((Safe_log10(y_23(ind_23, 16)) - Safe_log10(data_m(2, i))) ...
%             / max(Safe_log10(data_m(2, :)))) .^ 2;
%         err = err + ((Safe_log10(y_23(ind_23, 17)) - Safe_log10(data_nk(2, i))) ...
%             / max(Safe_log10(data_nk(2, :)))) .^ 2;
%         % err of ZT11
%         ind_11 = t == (data_m(1, i) + t_IAV_11);
%         err = err + ((Safe_log10(y_11(ind_11, 16)) - Safe_log10(data_m(3, i))) ...
%             / max(Safe_log10(data_m(3, :)))) .^ 2;
%         err = err + ((Safe_log10(y_11(ind_11, 17)) - Safe_log10(data_nk(3, i))) ...
%             / max(Safe_log10(data_nk(3, :)))) .^ 2;
%     end
%     
%     % cost function of T and TE
%     for i = 1:6
%         % err of ZT23
%         ind_23 = t == (data_te(1, i) + t_IAV_23);
%         err = err + ((Safe_log10(y_23(ind_23, 19)) - Safe_log10(data_te(2, i))) ...
%             / max(Safe_log10(data_te(2, :)))) .^ 2;
%         % err of ZT11
%         ind_11 = t == (data_te(1, i) + t_IAV_11);
%         err = err + ((Safe_log10(y_11(ind_11, 19)) - Safe_log10(data_te(3, i))) ...
%             / max(Safe_log10(data_te(3, :)))) .^ 2;
%     end

    % cost function of M and NK - non log error
    for i = 1:4
        % err of ZT23
        ind_23 = t == (data_m(1, i) + t_IAV_23);
        err = err + ((y_23(ind_23, 16) - data_m(2, i)) ...
            / max(data_m(2, :))) .^ 2;
        err = err + ((y_23(ind_23, 17) - data_nk(2, i)) ...
            / max(data_nk(2, :))) .^ 2;
        % err of ZT11
        ind_11 = t == (data_m(1, i) + t_IAV_11);
        err = err + ((y_11(ind_11, 16) - data_m(3, i)) ...
            / max(data_m(3, :))) .^ 2;
        err = err + ((y_11(ind_11, 17) - data_nk(3, i)) ...
            / max(data_nk(3, :))) .^ 2;
    end
    
    % cost function of T and TE - non log error
    for i = 1:6
        % err of ZT23
        ind_23 = t == (data_te(1, i) + t_IAV_23);
        err = err + ((y_23(ind_23, 19) - data_te(2, i)) ...
            / max(data_te(2, :))) .^ 2;
        % err of ZT11
        ind_11 = t == (data_te(1, i) + t_IAV_11);
        err = err + ((y_11(ind_11, 19) - data_te(3, i)) ...
            / max(data_te(3, :))) .^ 2;
    end
    
end