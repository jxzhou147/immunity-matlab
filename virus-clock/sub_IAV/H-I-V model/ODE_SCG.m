function dydt = ODE_SCG(t, y, par_scg)

    par_cell = num2cell(par_scg);
    [beta, d_I, p_V, c_V] = deal(par_cell{:});

    y_cell = num2cell(y);
    [H, I, V] = deal(y_cell{:});

    dydt = zeros(3, 1);     % H, I_1, I_2, V

    % ODEs
    dydt(1) = -1 * beta * H * V;

    dydt(2) = beta * H * V - d_I * I;

    dydt(3) = p_V * I - c_V * V;

end