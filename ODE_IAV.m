function dydt = ODE_IAV(t, y, par)

    dydt = zeros(8, 1);
    
    % parameters
    par = num2cell(par);
    [beta, gamma, eta, alpha_MH, alpha_MI, alpha_KI, alpha_TI, ...
        alpha_MV, c_IM, c_DM, c_M, c_IK, c_MK, c_MT, n_IM, n_DM, ...
        n_IK, n_MK, n_MT, K_IM, K_DM, K_M, K_IK, K_MK, K_MT, K_T, ...
        d_I, d_M, d_K, d_T] = deal(par{:});
    
    % variables
    y_cell = num2cell(y);
    [H, If, V, D, M, K, T, T_E] = deal(y_cell{:});
    
    % equations
    % dH/dt
    if H < 0
        H = 0;
        dydt(1) = 0;
    else
        dydt(1) = -1 * beta * H * V - alpha_MH * M * H;
    end

    % dI/dt
    if If < 0
        If = 0;
        dydt(2) = 0;
    else
        dydt(2) = beta * H * V - alpha_MI * M * If - ...
            alpha_KI * K * If - alpha_TI * T_E * If - d_I * If;
    end

    if V < 0
        V = 0;
        dydt(3) = 0;
    else
        dydt(3) = gamma * If - alpha_MV * M * V;
    end

    if D < 0
        D = 0;
        dydt(4) = 0;
    else
        dydt(4) = alpha_MH * M * H + alpha_MI * M * If + ...
            alpha_KI * K * If + alpha_TI * T_E * If + d_I * If;
    end

    if M < 0
        M = 0;
        dydt(5) = 0;
    else
        dydt(5) = c_IM * FracNoInf(RealRootPromise(If, n_IM), (RealRootPromise(K_IM, n_IM) + RealRootPromise(If, n_IM))) ... 
            + c_DM * FracNoInf(RealRootPromise(D, n_DM), (RealRootPromise(K_DM, n_DM) + RealRootPromise(D, n_DM))) ...
            + c_M * FracNoInf(M, (K_M + M)) - d_M * M;
    end

    if K < 0
        K = 0;
        dydt(6) = 0;
    else
        dydt(6) = c_IK * FracNoInf(RealRootPromise(If, n_IK), (RealRootPromise(K_IK, n_IK) + RealRootPromise(If, n_IK))) + ...
            c_MK * FracNoInf(RealRootPromise(M, n_MK), (RealRootPromise(K_MK, n_MK) + RealRootPromise(M, n_MK))) - d_K * K;
    end

    if T < 0
        T = 0;
        dydt(7) = 0;
    else
        dydt(7) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
            FracNoInf(RealRootPromise(M, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(M, n_MT))) * T;
    end

    if T_E < 0
        T_E = 0;
        dydt(8) = 0;
    else
        dydt(8) = c_MT * FracNoInf(RealRootPromise(M, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(M, n_MT))) * T - d_T * T_E;
    end
    
end