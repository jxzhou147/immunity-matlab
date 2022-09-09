function dydt = adjusted_ODE_IAV(t, y, par_IAV)
    
    dydt = zeros(15, 1);
    
    % variables of IAV model
    y_cell = num2cell(y);
    [H, If, D, V, M_0, M, Mono, N, IL1b, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
    % parameters of IAV model
    par_IAV = num2cell(par_IAV);
    [beta, gamma, eta, ...
        a_NH, a_MH, a_MI, a_NI, a_KI, a_TI, a_MD, a_NV, ...
        c_IL1b_M, c_IL1b_Mono, c_IL1b_N, c_IL1b_K, c_IM, c_CCL2_Mono, c_N, c_CXCL5_N, c_M_IL1b, c_I_IL1b, c_D_IL1b, ...
        c_M_IL10, c_M_CCL2, c_I_CCL2, c_N_CXCL5, c_I_CXCL5, c_K, c_MK, c_IK, c_MT, ...
        K_IL1b_M, K_IL1b_Mono, K_IL1b_N, K_IL1b_K, K_MI, K_I_M, K_D_IL1b, K_IL10_IL1b, K_IL10_CCL2, K_IL10_CXCL5, K_T, K_MT, n_I_M, n_D_IL1b, n_MT, ...
        d_H, d_I, d_V, d_M0, d_M, d_Mono, d_N, d_IL1b, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, ...
        Mono_tot, b_H, b_M0, b_N, b_K] = deal(par_IAV{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % equations of clock-controlled IAV model
    dydt(1) = (a_NH * b_N + a_MH * b_M0 + d_H) * b_H - beta * H * V - ...
        a_NH * N * H - a_MH * (M + Mono) * H - d_H * H;
%     dydt(1) = d_H * b_H - beta * H * V - a_NH * (N - b_N) * H - d_H * H;
%     dydt(1) = - beta * H * V;

    dydt(2) = beta * H * V - a_MI * FracNoInf((M + Mono), (K_MI + M + Mono)) * If - ...
        a_NI * N * If - a_KI * K * If - a_TI * T_E * If - d_I * If;
   
    dydt(3) = a_NI * N * If + a_KI * K * If + a_TI * T_E * If + d_I * If + ...
        a_NH * N * H + a_MH * (M + Mono) * H - a_MD * (M + Mono) * D;

%     if (t > 15)
%         gamma = (1 - exp(-0.05 * (t - 15))) * gamma;
%     end
    dydt(4) = gamma * If - a_NV * N * V - d_V * V;

    dydt(5) = b_M0 * d_M0 - (1 + c_IL1b_M * FracNoInf(IL1b, K_IL1b_M +IL1b)) * ...
        c_IM * FracNoInf(RealRootPromise(If, n_I_M), RealRootPromise(K_I_M, n_I_M) + RealRootPromise(If, n_I_M)) * M_0 - ...
        d_M0 * M_0;
    
    dydt(6) = (1 + c_IL1b_M * FracNoInf(IL1b, K_IL1b_M +IL1b)) * ...
        c_IM * FracNoInf(RealRootPromise(If, n_I_M), RealRootPromise(K_I_M, n_I_M) + RealRootPromise(If, n_I_M)) * M_0 - ...
        d_M * M;
    
%     dydt(7) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * CCL2 * (Mono_tot - Mono) - d_Mono * Mono;
    dydt(7) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * 1000 * FracNoInf(CCL2, CCL2 + 500) * (Mono_tot - Mono) - d_Mono * Mono;

    dydt(8) = (c_N + c_CXCL5_N * (1 + c_IL1b_N * FracNoInf(IL1b, (K_IL1b_N + IL1b))) * CXCL5) * ((c_N + d_N) / c_N * b_N - N) - d_N * N;

    if (t > 15)
        c_M_IL1b = (1 - exp(-0.05 * (t - 15))) * c_M_IL1b;
        c_I_IL1b = (1 - exp(-0.05 * (t - 15))) * c_I_IL1b;
    end
    dydt(9) = (FracNoInf(K_IL10_IL1b, K_IL10_IL1b + IL10) * c_M_IL1b * (M + Mono) + c_I_IL1b * If) * ...
        (1 + c_D_IL1b * FracNoInf(RealRootPromise(D, n_D_IL1b), RealRootPromise(K_D_IL1b, n_D_IL1b) + RealRootPromise(D, n_D_IL1b))) - d_IL1b * IL1b;

    dydt(10) = c_M_IL10 * (M + Mono) - d_IL10 * IL10;

    if (t > 15)
        c_M_CCL2 = (1 - exp(-0.05 * (t - 15))) * c_M_CCL2;
        c_I_CCL2 = (1 - exp(-0.05 * (t - 15))) * c_I_CCL2;
    end
    dydt(11) = FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * (M + Mono) + ...
        c_I_CCL2 * If - d_CCL2 * CCL2;

    if (t > 15)
        c_N_CXCL5 = (1 - exp(-0.05 * (t - 15))) * c_N_CXCL5;
        c_I_CXCL5 = (1 - exp(-0.05 * (t - 15))) * c_I_CXCL5;
    end
    dydt(12) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_N_CXCL5 * N + ...
        c_I_CXCL5 * If - d_CXCL5 * CXCL5;
    
    dydt(13) = (1 + c_IL1b_K * FracNoInf(IL1b, (K_IL1b_K + IL1b))) * (c_K + c_MK * M + c_IK * If) * ((c_K + d_K) / c_K * b_K - K) - d_K * K;

%    dydt(13) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
%         FracNoInf(RealRootPromise((M + Mono), n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise((M + Mono), n_MT))) * T;
% 
%     dydt(14) = c_MT * FracNoInf(RealRootPromise((M + Mono), n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise((M + Mono), n_MT))) * T - d_T * T_E;
    dydt(14) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
        FracNoInf(RealRootPromise(K, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(K, n_MT))) * T;

    dydt(15) = c_MT * FracNoInf(RealRootPromise(K, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(K, n_MT))) * T - d_T * T_E;

    end
