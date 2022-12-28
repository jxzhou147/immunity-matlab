function dydt = ODE_SCG(t, y, par_scg)

    par_cell = num2cell(par_scg);
    [beta, tau_E, tau_I, n_E, n_I, p_inf, p_tot, c_inf, c_tot, ...
        p_IL6, p_IL1b, p_IL10, d_IL6, d_IL1b, d_IL10] = deal(par_cell{:});

    dydt = zeros(n_E + n_I + 6, 1);     % T, Ei, Ij, V_inf, V_tot, IL6, IL1b, IL10

    % T
    dydt(1) = -1 * beta * y(1) * y(n_E+n_I+2);

    % ODEs for E_1,...,E_{n_E}
    dydt(2) = beta * y(1) * y(n_E+n_I+2) - 1 * n_E / tau_E * y(2);
    for i = 3:(n_E+1)
        dydt(i) = n_E / tau_E * y(i-1) - n_E / tau_E * y(i);
    end

    % ODEs for I_1,...,I_{n_I}
    dydt(n_E+2) = n_E / tau_E * y(n_E+1) - n_I / tau_I * y(n_E+2);
    for j = (n_E+3):(n_E+n_I+1)
        dydt(j) = n_I / tau_I * y(j-1) - n_I / tau_I * y(j);
    end

    % V_inf, V_tot, IL6, IL1b, IL10
    for j = (n_E+2):(n_E+n_I+1)
        dydt(n_E+n_I+2) = dydt(n_E+n_I+2) + p_inf * y(j);
    end
    dydt(n_E+n_I+2) = dydt(n_E+n_I+2) - c_inf * y(n_E+n_I+2);

    for j = (n_E+2):(n_E+n_I+1)
        dydt(n_E+n_I+3) = dydt(n_E+n_I+3) + p_tot * y(j);
    end
    dydt(n_E+n_I+3) = dydt(n_E+n_I+3) - c_tot * y(n_E+n_I+3);

    for j = (n_E+2):(n_E+n_I+1)
        dydt(n_E+n_I+4) = dydt(n_E+n_I+4) + p_IL6 * y(j);
    end
    dydt(n_E+n_I+4) = dydt(n_E+n_I+4) - d_IL6 * y(n_E+n_I+4);

    for j = (n_E+2):(n_E+n_I+1)
        dydt(n_E+n_I+5) = dydt(n_E+n_I+5) + p_IL1b * y(j);
    end
    dydt(n_E+n_I+5) = dydt(n_E+n_I+5) - d_IL1b * y(n_E+n_I+5);

    for j = (n_E+2):(n_E+n_I+1)
        dydt(n_E+n_I+6) = dydt(n_E+n_I+6) + p_IL10 * y(j);
    end
    dydt(n_E+n_I+6) = dydt(n_E+n_I+6) - d_IL10 * y(n_E+n_I+6);

end