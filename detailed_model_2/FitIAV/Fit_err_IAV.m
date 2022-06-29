function err = Fit_err_IAV(data_h, data_m, data_mono, data_neu, data_nk, data_v, data_t, data_te, data_il1b, data_ccl2, parFit, parBase)
    
    % meaning of inputs
    % data*(1): time; data*(2): value of ZT23; data*(3): value of ZT11
    % parFit: 74 parameters
    
    % fix some parameters
%     parFit(1) = -6.633766;
%     parFit(2) = 3.247969;
%     parFit(3) = -0.124729;
    
    % translate parFit to par in the odes
    log_par_ind = [1:43 47:60];
    par_IAV = [parFit; parBase];
    for i = log_par_ind
        par_IAV(i) = 10 .^ parFit(i);
    end
         
    % solve odes
    tmax = 250;
    tspan = 0:1:tmax;
    y0 = zeros(15, 1);
    y0(4) = 40;
    y0(1) = parBase(1);
    y0(5) = parBase(2);
    y0(8) = parBase(3);
    y0(13) = parBase(4);
    y0(14) = 10;
    
    
    % For ZT23
%     [t, y_23] = ode15s(@ODE_IAV, tspan, y0, [], par_IAV);
%     y_23 = real(y_23);
    
    % For ZT11
    [t, y_11] = ode15s(@ODE_IAV, tspan, y0, [], par_IAV);
    y_11 = real(y_11);

    err = 0;
    
    data_ind = 3;   % 2: ZT23, 3: ZT11
    
    % cost function of V - log error
    for i = 1:8
        % err of ZT23
%         ind_23 = t == data_v(1, i);
%         err = err + ((y_23(ind_23, 4) - data_v(2, i)) ...
%             / max(data_v(2, :))) .^ 2;
        % err of ZT11
        ind_11 = t == data_v(1, i);
        err = err + ((y_11(ind_11, 4) - data_v(data_ind, i)) ...
            / max(data_v(data_ind, :))) .^ 2;
    end
    

    % cost function of H, M, Neu and NK - non log error
%     for i = 1:4
%         % err of ZT23
%         ind_23 = t == data_m(1, i);
%         err = err + ((y_23(ind_23, 1) - data_h(2, i)) ...
%             / max(data_h(2, :))) .^ 2;
%         err = err + ((y_23(ind_23, 5) - data_m(2, i)) ...
%             / max(data_m(2, :))) .^ 2;
%         err = err + ((y_23(ind_23, 6) - data_neu(2, i)) ...
%             / max(data_neu(2, :))) .^ 2;
%         err = err + ((y_23(ind_23, 11) - data_nk(2, i)) ...
%             / max(data_nk(2, :))) .^ 2;
%     end
    
    for i = 1:3
        % err of ZT11
        ind_11 = t == data_m(1, i);
        err = err + ((y_11(ind_11, 1) - data_h(data_ind, i)) ...
            / max(data_h(data_ind, :))) .^ 2;
        err = err + ((y_11(ind_11, 5) + y_11(ind_11, 6) - data_m(data_ind, i)) ...
            / max(data_m(data_ind, :))) .^ 2;
        err = err + ((y_11(ind_11, 7) - data_mono(data_ind, i)) ...
            / max(data_mono(data_ind, :))) .^ 2;
        err = err + ((y_11(ind_11, 8) - data_neu(data_ind, i)) ...
            / max(data_neu(data_ind, :))) .^ 2;
        err = err + ((y_11(ind_11, 13) - data_nk(data_ind, i)) ...
            / max(data_nk(data_ind, :))) .^ 2;
    end
    
    % cost function of T and TE - non log error
    for i = 1:6
        % err of ZT23
%         ind_23 = t == data_t(1, i);
% %         err = err + ((y_23(ind_23, 12) - data_t(2, i)) ...
% %             / max(data_t(2, :))) .^ 2;
%         err = err + ((y_23(ind_23, 13) - data_te(2, i)) ...
%             / max(data_te(2, :))) .^ 2;
        
        % err of ZT11
        ind_11 = t == data_t(1, i);
%         err = err + ((y_23(ind_23, 12) - data_t(2, i)) ...
%             / max(data_t(2, :))) .^ 2;
        err = err + ((y_11(ind_11, 15) - data_te(data_ind, i)) ...
            / max(data_te(data_ind, :))) .^ 2;
    end
    
    % cost function of IL6 and CCL2 - non log error
    for i = 1:4
%         % err of ZT23
%         ind_23 = t == data_il1b(1, i);
%         err = err + ((y_23(ind_23, 7) - data_il1b(2, i)) ...
%             / max(data_il1b(2, :))) .^ 2;
%         err = err + ((y_23(ind_23, 9) - data_ccl2(2, i)) ...
%             / max(data_ccl2(2, :))) .^ 2;
        
        % err of ZT11
        ind_11 = t == data_ccl2(1, i);
%         err = err + ((y_11(ind_11, 9) - data_il1b(data_ind, i)) ...
%             / max(data_il1b(data_ind, :))) .^ 2;
        err = err + ((y_11(ind_11, 11) - data_ccl2(data_ind, i)) ...
            / max(data_ccl2(data_ind, :))) .^ 2;
    end

end