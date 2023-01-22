function dydt = ODE_IAV(t, y, par_IAV)
    
    dydt = zeros(13, 1);
    
    % variables of IAV model
    y_cell = num2cell(y);
    [H, I1, I2, V, Mono, N, IL1b, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
    % parameters of IAV model
    par_IAV = num2cell(par_IAV);
    [zeta, beta, k_I, gamma, eta, r, a_MI, a_NI, a_KI, a_TI, ...
        c_IL1b_Mono, c_IL1b_N, c_CCL2_Mono, c_CXCL5_N, c_M_IL1b, c_I_IL1b, ...
        c_M_IL10, c_M_CCL2, c_I_CCL2, c_N_CXCL5, c_I_CXCL5, c_I_K, c_M_T, ...
        K_V, K_IL1b_Mono, K_IL1b_N, K_CCL2_Mono, K_CXCL5_N, K_IL10_IL1b, K_IL10_CCL2, K_IL10_CXCL5, ...
        K_I_K, K_T, K_M_T, n_CCL2_Mono, n_CXCL5_N, n_I_K, n_M_T, ...
        d_I, d_V, d_Vs, d_Mono, d_N, d_IL1b, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, b_H, b_K] ...
        = deal(par_IAV{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dydt(1) = zeta * H * (1 - H / b_H) - beta * H * V;

    dydt(2) = beta * H * V - k_I * I1;
    
    dydt(3) = k_I * I1 - a_MI * Mono * I2 - a_NI * N * I2 - a_KI  * K * I2 - a_TI * T_E * I2 - d_I * I2;

%     dydt(4) = gamma * I2 - d_V * V - d_Vs * V / (K_V + V);
%     dydt(4) = gamma * I2 - T_E / (T_E + 150) *  d_V * V - d_Vs * V / (K_V + V);
    if (t <= 1)
        dydt(4) = gamma * I2 - T_E / (T_E + 150) *  d_V * V - d_Vs * V / (K_V + V) + 100;
    else
        dydt(4) = gamma * I2 - T_E / (T_E + 150) *  d_V * V - d_Vs * V / (K_V + V);
    end
    dydt(5) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * ...
        FracNoInf(RealRootPromise(CCL2, n_CCL2_Mono), RealRootPromise(K_CCL2_Mono, n_CCL2_Mono) + RealRootPromise(CCL2, n_CCL2_Mono)) - d_Mono * Mono;

    dydt(6) = c_CXCL5_N * (1 + c_IL1b_N * FracNoInf(IL1b, (K_IL1b_N + IL1b))) * ...
        FracNoInf(RealRootPromise(CXCL5, n_CXCL5_N), RealRootPromise(K_CXCL5_N, n_CXCL5_N) + RealRootPromise(CXCL5, n_CXCL5_N)) - d_N * N;

    dydt(7) = FracNoInf(K_IL10_IL1b, K_IL10_IL1b + IL10) * c_M_IL1b * Mono + ...
        c_I_IL1b * I2 - d_IL1b * IL1b;

    dydt(8) = c_M_IL10 * Mono - d_IL10 * IL10;

    dydt(9) = FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * Mono + ...
        c_I_CCL2 * I2 - d_CCL2 * CCL2;

    dydt(10) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_N_CXCL5 * N + ...
        c_I_CXCL5 * I2 - d_CXCL5 * CXCL5;

    dydt(11) = c_I_K * FracNoInf(RealRootPromise(I2, n_I_K), RealRootPromise(K_I_K, n_I_K) + RealRootPromise(I2, n_I_K)) - d_K * (K - b_K);

    dydt(12) = eta * T * (1 - FracNoInf(T, K_T)) - c_M_T * ...
    FracNoInf(RealRootPromise(Mono, n_M_T), (RealRootPromise(K_M_T, n_M_T) + RealRootPromise(Mono, n_M_T))) * T;

    dydt(13) = r * c_M_T * FracNoInf(RealRootPromise(Mono, n_M_T), (RealRootPromise(K_M_T, n_M_T) + RealRootPromise(Mono, n_M_T))) * T - d_T * T_E;

end