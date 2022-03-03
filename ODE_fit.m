function err = ODE_fit(data_v, data_m, data_nk, data_t, data_te, parFit)

    % translate parFit to par in the odes
    par = parFit;
    for i = [(1:14), (20:30)]
        par(i) = 10 .^ parFit(i);
    end
    
    % solve odes
    tmax = 250;
    tspan = 0:1:tmax;
    y0 = [10^7, 40, 0, 0, 0, 0, 0, 0];
    [t, y] = ode45(@ODE_IAV, tspan, y0, [], par);

    err = 0;
    
    % cost function of V
    for i = 1:8
        ind = t == data_v(1, i);
        err = err + ((Safe_log10(y(ind, 3)) - Safe_log10(data_v(2, i))) ...
            / max(Safe_log10(data_v(2, :)))) ^ 2;
    end
    
    % cost function of M and NK
    for i = 1:4
        ind = t == data_m(1, i);
        err = err + ((Safe_log10(y(ind, 5)) - Safe_log10(data_m(2, i))) ...
            / max(Safe_log10(data_m(2, :)))) ^ 2;
        err = err + ((Safe_log10(y(ind, 6)) - Safe_log10(data_nk(2, i))) ...
            / max(Safe_log10(data_nk(2, :)))) ^2;
    end
     
    % cost function of T and TE
%     for i = 1:2
%         ind = t == data_t(1, i);
%         err = err + ((y(ind, 7) - data_t(2, i)) / max(data_t(2, :))) ^ 2;
%         err = err + ((y(ind, 8) - data_te(2, i)) / max(data_te(2, :))) ^ 2;
%     end

end