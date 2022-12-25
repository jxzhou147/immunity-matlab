function dydt = ODE_SCG(t, y, par_scg)

    par_cell = num2cell(par_scg);
    [tau_E, tau_I, n_E, n_I, p_inf, p_tot, c_inf, c_B] = deal(par_cell{:});

    dydt = zeros(n_E + n_I + 2, 1);     % Ei, Ij, V

    % ODEs for E_1,...,E_{n_E}
    dydt(1) = -1 * n_E / tau_E * y(1);
    for i = 2:n_E
        dydt(i) = n_E / tau_E * y(i-1) - n_E / tau_E * y(i);
    end

    % ODEs for I_1,...,I_{n_I}
    dydt(n_E+1) = n_E / tau_E * y(n_E) - n_I / tau_I * y(n_E+1);
    for j = (n_E+2):(n_E+n_I)
        dydt(j) = n_I / tau_I * y(j-1) - n_I / tau_I * y(j);
    end

    for j = (n_E+1):(n_E+n_I)
        dydt(n_E+n_I+1) = dydt(n_E+n_I+1) + p_inf * y(j);
    end
    dydt(n_E+n_I+1) = dydt(n_E+n_I+1) - c_inf * y(n_E+n_I+1);

    for j = (n_E+1):(n_E+n_I)
        dydt(n_E+n_I+2) = dydt(n_E+n_I+2) + p_tot * y(j);
    end
    dydt(n_E+n_I+2) = dydt(n_E+n_I+2) - c_B * y(n_E+n_I+2);

end