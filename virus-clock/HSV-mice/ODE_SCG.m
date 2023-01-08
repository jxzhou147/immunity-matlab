function dydt = ODE_SCG(t, y, par_scg)

    par_cell = num2cell(par_scg);
    [beta, k_I, d_I, p_V, c_V] = deal(par_cell{:});

    y_cell = num2cell(y);
    [H, I1, I2, V] = deal(y_cell{:});

    dydt = zeros(4, 1);     % H, I_1, I_2, V

    % ODEs
    dydt(1) = -1 * beta * H * V;

    dydt(2) = beta * H * V - k_I * I1;

    dydt(3) = k_I * I1 - d_I * I2;

    dydt(4) = p_V * I2 - c_V * V;

end