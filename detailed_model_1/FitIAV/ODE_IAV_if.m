function dydt = ODE_IAV(t, y, par_IAV)
    
    dydt = zeros(13, 1);
    
    % variables of IAV model
    y_cell = num2cell(y);
    [H, If, D, V, M, N, IL6, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
    % parameters of IAV model
    par_IAV = num2cell(par_IAV);
    [beta, gamma, eta, ...
        a_NH, a_MI, a_NI, a_KI, a_TI, a_NV, ...
        c_IL6_M, c_CCL2_M, c_CXCL5_N, c_M_IL6, c_K_IL6, c_N_IL6, c_I_IL6, ...
        c_IM, c_DM, c_M_IL10, c_M_CCL2, c_I_CCL2, c_M_CXCL5, c_I_CXCL5, c_IK, c_CCL2_K, c_MT, ...
        K_IL6_M, K_IM, K_DM, K_CCL2_M, K_CXCL5_N, K_IL10_IL6, K_IL10_CCL2, ...
        K_IL10_CXCL5, K_IK, K_CCL2_K, K_T, K_MT, ...
        n_IM, n_DM, n_CCL2_M, n_CXCL5_N, n_IK, n_CCL2_K, n_MT, ...
        d_H, d_I, d_V, d_M, d_N, d_IL6, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, ...
        b_H, b_M, b_N, b_K] = deal(par_IAV{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % equations of clock-controlled IAV model
    if (H < 0)
        H = 0;
        dydt(1) = 0;
    else
        dydt(1) = d_H * b_H - beta * H * V - a_NH * (N - b_N) * H - d_H * H;
%         dydt(1) = -beta * H * V;
    end

    if (If < 0)
        If = 0;
        dydt(2) = 0;
    else
        dydt(2) = beta * H * V - (1 + c_IL6_M * FracNoInf(IL6, (K_IL6_M + IL6))) * a_MI * M * If - ...
            a_NI * N * If - a_KI * K * If - a_TI * T_E * If - d_I * If;
    end

    if (D < 0)
        D = 0;
        dydt(3) = 0;
    else
        dydt(3) = (1 + c_IL6_M * FracNoInf(IL6, (K_IL6_M + IL6))) * a_MI * M * If + ...
            a_NI * N * If + a_KI * K * If + a_TI * T_E * If + d_I * If + ...
            a_NH * (N - b_N) * H + d_H * H;
%         dydt(3) = (1 + c_IL6_M * FracNoInf(IL6, (K_IL6_M + IL6))) * a_MI * M * If + ...
%             a_NI * N * If + a_KI * K * If + a_TI * T_E * If + d_I * If;
    end

    if (V < 0)
        V = 0;
        dydt(4) = 0;
    else
        dydt(4) = gamma * If - a_NV * N * V - d_V * V;
    end

    if (M < 0)
        M = 0;
        dydt(5) = 0;
    else
        dydt(5) = b_M * d_M + c_IM * FracNoInf(RealRootPromise(If, n_IM), (RealRootPromise(K_IM, n_IM) + RealRootPromise(If, n_IM))) + ...
            c_DM * FracNoInf(RealRootPromise(D, n_DM), (RealRootPromise(K_DM, n_DM) + RealRootPromise(D, n_DM))) + ...
            c_CCL2_M * FracNoInf(RealRootPromise(CCL2, n_CCL2_M), (RealRootPromise(K_CCL2_M, n_CCL2_M) + RealRootPromise(CCL2, n_CCL2_M))) - ...
            d_M * M;
    end

    if (N < 0)
        N = 0;
        dydt(6) = 0;
    else
        dydt(6) = b_N * d_N + c_CXCL5_N * FracNoInf(RealRootPromise(CXCL5, n_CXCL5_N), (RealRootPromise(K_CXCL5_N, n_CXCL5_N) + RealRootPromise(CXCL5, n_CXCL5_N))) - ...
            d_N * N;
    end

    if (IL6 < 0)
        IL6 = 0;
        dydt(7) = 0;
    else
        dydt(7) = FracNoInf(K_IL10_IL6, K_IL10_IL6 + IL10) * c_M_IL6 * M + ...
            c_K_IL6 * K + c_N_IL6 * N + c_I_IL6 * If - d_IL6 * IL6;
    end

    if (IL10 < 0)
        IL10 = 0;
        dydt(8) = 0;
    else
        dydt(8) = c_M_IL10 * M - d_IL10 * IL10;
    end

    if (CCL2 < 0)
        CCL2 = 0;
        dydt(9) = 0;
    else
        dydt(9) = FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * M + ...
            c_I_CCL2 * If - d_CCL2 * CCL2;
    end

    if (CXCL5 < 0)
        CXCL5 = 0;
        dydt(10) = 0;
    else
        dydt(10) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_M_CXCL5 * M + ...
            c_I_CXCL5 * If - d_CXCL5 * CXCL5;
    end

    if (K < 0)
        K = 0;
        dydt(11) = 0;
    else
        dydt(11) = b_K * d_K + c_IK * FracNoInf(RealRootPromise(If, n_IK), (RealRootPromise(K_IK, n_IK) + RealRootPromise(If, n_IK))) + ...
            c_CCL2_K * FracNoInf(RealRootPromise(CCL2, n_CCL2_K), RealRootPromise(K_CCL2_K, n_CCL2_K) + RealRootPromise(CCL2, n_CCL2_K)) ...
            - d_K * K;
    end

    if (T < 0)
        T = 0;
        dydt(12) = 0;
    else
        dydt(12) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
            FracNoInf(RealRootPromise(M, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(M, n_MT))) * T;
    end

    if (T_E < 0)
        T_E = 0;
        dydt(13) = 0;
    else
        dydt(13) = c_MT * FracNoInf(RealRootPromise(M, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(M, n_MT))) * T - d_T * T_E;
    end
    
end