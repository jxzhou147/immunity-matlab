function dydt = ODE_IAV(t, y, par_IAV)

    dydt = zeros(14, 1);    % y(16):weight change
    
    % variables of IAV model
    y_cell = num2cell(y);
    [H, If, V, M_0, M, Mono, N, IL1b, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
    % parameters of IAV model
    par_IAV = num2cell(par_IAV);
    [zeta, beta, gamma, eta, r, ...
        a_MI, a_NI, a_KI, a_TI, a_NV, ...
        c_IL1b_M, c_IL1b_Mono, c_IL1b_N, c_IM, c_CCL2_Mono, c_CXCL5_N, c_M_IL1b, c_I_IL1b, ...
        c_M_IL10, c_M_CCL2, c_I_CCL2, c_N_CXCL5, c_I_CXCL5, c_IK, c_MT, ...
        K_V, K_IL1b_M, K_IL1b_Mono, K_IL1b_N, K_MI, K_I_M, K_CCL2_Mono, K_CXCL5_N, K_IL10_IL1b, K_IL10_CCL2, K_IL10_CXCL5, K_IK, K_T, K_MT, n_I_M, n_CCL2_Mono, n_CXCL5_N, n_IK, n_MT, ...
        d_I, d_V, d_Vs, d_M0, d_M, d_Mono, d_N, d_IL1b, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, ...
        b_H, b_M0, b_N, b_K] = deal(par_IAV{:});
    
    % parameters of inflammation and weigth change
%     par_infla = num2cell(par_infla);
%     [sigma_M, sigma_Mono, sigma_N, sigma_I, sigma_Infla, K_I, K_In, n_I, n_In, d_wc] = deal(par_infla{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % equations of IAV model

    dydt(1) = zeta * H * (1 - H / b_H ) - beta * H * V;
    
    dydt(2) = beta * H * V - a_MI * FracNoInf((M + Mono), (K_MI + M + Mono)) * If - ...
        a_NI * N * If - a_KI * K * If - a_TI * T_E * If - d_I * If;

%     dydt(3) = gamma * If - a_NV * N * V - d_V * V - 8e-2 * V / (1e-2 + V);
    dydt(3) = gamma * If - a_NV * N * V - d_V * V - d_Vs * V / (K_V + V);

    dydt(4) = b_M0 * d_M0 - (1 + c_IL1b_M * FracNoInf(IL1b, K_IL1b_M +IL1b)) * ...
        c_IM * FracNoInf(RealRootPromise(If, n_I_M), RealRootPromise(K_I_M, n_I_M) + RealRootPromise(If, n_I_M)) * M_0 - ...
        d_M0 * M_0;

    dydt(5) = (1 + c_IL1b_M * FracNoInf(IL1b, K_IL1b_M +IL1b)) * ...
        c_IM * FracNoInf(RealRootPromise(If, n_I_M), RealRootPromise(K_I_M, n_I_M) + RealRootPromise(If, n_I_M)) * M_0 - ...
        d_M * M;

    dydt(6) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * ...
        FracNoInf(RealRootPromise(CCL2, n_CCL2_Mono), RealRootPromise(K_CCL2_Mono, n_CCL2_Mono) + RealRootPromise(CCL2, n_CCL2_Mono)) - d_Mono * Mono;
    
    dydt(7) = c_CXCL5_N * (1 + c_IL1b_N * FracNoInf(IL1b, (K_IL1b_N + IL1b))) * ...
        FracNoInf(RealRootPromise(CXCL5, n_CXCL5_N), RealRootPromise(K_CXCL5_N, n_CXCL5_N) + RealRootPromise(CXCL5, n_CXCL5_N)) - d_N * (N - b_N);

    dydt(8) = (FracNoInf(K_IL10_IL1b, K_IL10_IL1b + IL10) * c_M_IL1b * (M + Mono) + c_I_IL1b * If) - d_IL1b * IL1b;

    dydt(9) = c_M_IL10 * (M + Mono) - d_IL10 * IL10;

    dydt(10) = FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * (M + Mono) + ...
        c_I_CCL2 * If - d_CCL2 * CCL2;

    dydt(11) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_N_CXCL5 * (N - b_N) + ...
        c_I_CXCL5 * If - d_CXCL5 * CXCL5;
    
    dydt(12) = c_IK * FracNoInf(RealRootPromise(If, n_IK), RealRootPromise(K_IK, n_IK) + RealRootPromise(If, n_IK)) - d_K * (K - b_K);

    dydt(13) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
        FracNoInf(RealRootPromise(Mono, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(Mono, n_MT))) * T;

    dydt(14) = r * c_MT * FracNoInf(RealRootPromise(Mono, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(Mono, n_MT))) * T - d_T * T_E;
    
    %%%%%%%%%%%%%%%% dynamics of inflammation and weigth change
%     infla = sigma_M * M + sigma_Mono * Mono + sigma_N * N;
%     infla = sigma_Mono * Mono + sigma_N * (N - b_N);
%     dydt(16) = sigma_I * FracNoInf(RealRootPromise(If, n_I), RealRootPromise(K_I, n_I) + RealRootPromise(If, n_I)) + ...
%             sigma_Infla * FracNoInf(RealRootPromise(infla, n_In), RealRootPromise(K_In, n_In) + RealRootPromise(infla, n_In)) - d_wc * wc;
    

end