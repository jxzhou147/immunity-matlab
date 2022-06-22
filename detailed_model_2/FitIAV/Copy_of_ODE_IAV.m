function dydt = ODE_IAV(t, y, par_IAV)
    
    dydt = zeros(14, 1);
    
    % variables of IAV model
    y_cell = num2cell(y);
    [H, If, D, V, M, Mono, N, IL6, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
    % parameters of IAV model
    par_IAV = num2cell(par_IAV);
    [beta, gamma, eta, ...
        a_MI, a_NI, a_KI, a_TI, a_MD, a_NV, ...
        c_IL6_M, c_IL6_Mono, c_IL6_N, c_IL6_K, c_IM, c_CCL2_Mono, c_N, c_CXCL5_N, c_M_IL6, c_K_IL6, c_N_IL6, c_I_IL6, ...
        c_M_IL10, c_M_CCL2, c_I_CCL2, c_N_CXCL5, c_I_CXCL5, c_K, c_MK, c_IK, c_MT, ...
        K_IL6_M, K_IL6_Mono, K_IL6_N, K_IL6_K, K_MI, K_IL10_IL6, K_IL10_CCL2, K_IL10_CXCL5, K_T, K_MT, n_MT, ...
        d_H, d_I, d_V, d_M, d_Mono, d_N, d_IL6, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, ...
        Mono_tot, b_H, b_M, b_N, b_K] = deal(par_IAV{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % equations of clock-controlled IAV model
    dydt(1) = d_H * b_H - beta * H * V - d_H * H;
%     dydt(1) = d_H * b_H - beta * H * V - a_NH * (N - b_N) * H - d_H * H;
%     dydt(1) = - beta * H * V;

    dydt(2) = beta * H * V - (1 + c_IL6_M * FracNoInf(IL6, (K_IL6_M + IL6))) * a_MI * FracNoInf((M + Mono), (K_MI + M + Mono)) * If - ...
        a_NI * N * If - a_KI * K * If - a_TI * T_E * If - d_I * If;
   
    dydt(3) = (1 + c_IL6_M * FracNoInf(IL6, (K_IL6_M + IL6))) * a_MI * FracNoInf((M + Mono), (K_MI + M + Mono)) * If + ...
        a_NI * N * If + a_KI * K * If + a_TI * T_E * If + d_I * If -a_MD * (M + Mono) * D;

    dydt(4) = gamma * If - a_NV * N * V - d_V * V;

    dydt(5) = b_M * d_M - c_IM * M * If - d_M * M;
    
%     c_CCL2_Mono = (1 - exp(-0.1 * t)) * c_CCL2_Mono;
    dydt(6) = c_CCL2_Mono * (1 + c_IL6_Mono * FracNoInf(IL6, (K_IL6_Mono + IL6))) * CCL2 * (Mono_tot - Mono) - d_Mono * Mono;
    
    dydt(7) = (c_N + c_CXCL5_N * (1 + c_IL6_N * FracNoInf(IL6, (K_IL6_N + IL6))) * CXCL5) * ((c_N + d_N) / c_N * b_N - N) - d_N * N;

    dydt(8) = FracNoInf(K_IL10_IL6, K_IL10_IL6 + IL10) * c_M_IL6 * (M + Mono) + ...
        c_K_IL6 * K + c_N_IL6 * N + c_I_IL6 * If - d_IL6 * IL6;

    dydt(9) = c_M_IL10 * (M + Mono) - d_IL10 * IL10;

%     c_M_CCL2 = (1 - exp(-0.1 * t)) * c_M_CCL2;
    dydt(10) = FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * (M + Mono) + ...
        c_I_CCL2 * If - d_CCL2 * CCL2;

    dydt(11) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_N_CXCL5 * N + ...
        c_I_CXCL5 * If - d_CXCL5 * CXCL5;
    
    dydt(12) = (1 + c_IL6_K * FracNoInf(IL6, (K_IL6_K + IL6))) * (c_K + c_MK * (M + Mono) + c_IK * If) * ((c_K + d_K) / c_K * b_K - K) - d_K * K;

%    dydt(13) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
%         FracNoInf(RealRootPromise((M + Mono), n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise((M + Mono), n_MT))) * T;
% 
%     dydt(14) = c_MT * FracNoInf(RealRootPromise((M + Mono), n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise((M + Mono), n_MT))) * T - d_T * T_E;
    dydt(13) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
        FracNoInf(RealRootPromise(K, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(K, n_MT))) * T;

    dydt(14) = c_MT * FracNoInf(RealRootPromise(K, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(K, n_MT))) * T - d_T * T_E;

end