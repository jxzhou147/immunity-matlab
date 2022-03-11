function dydt = ODE_Clock_IAV(t, y, par_clock, par_IAV, t_IAV)
    
    dydt = zeros(19, 1);
    
    % variables of clock_IAV model
    y_cell = num2cell(y);
    [Per, Cry, Rev, Ror, Bmal1, PER, CRY, REV, ROR, ...
        BMAL1, PER_CRY, CLOCK_BMAL1, ...
        H, If, V, M, K, T, T_E] = deal(y_cell{:});
    
    % parameters of clock model
    % mRNA and protein degradation rate constants (in h^-1)
    par_clock1 = num2cell(par_clock(1:12));
    [dm_per, dm_cry, dm_rev, dm_ror, dm_bmal, dp_per, ...
        dp_cry, dp_rev, dp_ror, dp_bmal, d_pc, d_cb] = deal(par_clock1{:});
    
    % maximal transcription rates (in nmol l^-1 h^-1)
    par_clock2 = num2cell(par_clock(13:17));
    [vmax_per, vmax_cry, vmax_rev, vmax_ror, vmax_bmal] = deal(par_clock2{:});
    
    % activation ratios (demensionless)
    par_clock3 = num2cell(par_clock(18:22));
    [fold_per, fold_cry, fold_rev, fold_ror, fold_bmal] = deal(par_clock3{:});
    
    % regulation thresholds (in nmol/l)
    par_clock4 = num2cell(par_clock(23:33));
    [Ka_per_cb, Ki_per_pc, Ka_cry_cb, Ki_cry_pc, Ki_cry_rev, Ka_rev_cb, ...
        Ki_rev_pc, Ka_ror_cb, Ki_ror_pc, Ka_bmal_ror, Ki_bmal_rev] = deal(par_clock4{:});
    
    % hill coefficients (dimensionless)
    par_clock5 = num2cell(par_clock(34:44));
    [hill_per_cb, hill_per_pc, hill_cry_cb, hill_cry_pc, ...
        hill_cry_rev, hill_rev_cb, hill_rev_pc, hill_ror_cb, ...
        hill_ror_pc, hill_bmal_ror, hill_bmal_rev] = deal(par_clock5{:});
    
    % translation rates (in molecules per hour per mRNA)
    par_clock6 = num2cell(par_clock(45:49));
    [kp_per, kp_cry, kp_rev, kp_ror, kp_bmal] = deal(par_clock6{:});
    
    % complexation kinetic rates
    par_clock7 = num2cell(par_clock(50:53));
    [kass_cb, kass_pc, kdiss_cb, kdiss_pc] = deal(par_clock7{:});
    
    
    % parameters of IAV model
    par_IAV = num2cell(par_IAV);
    [beta, gamma, eta, alpha_MI, alpha_KI, alpha_TI, ...
        alpha_MV, c_IM, c_M, c_IK, c_MK, c_MT, n_IM, ...
        n_IK, n_MK, n_MT, K_IM, K_M, K_IK, K_MK, K_MT, K_T, ...
        K_TI, d_I, d_M, d_K, d_T, K_B, K_R] = deal(par_IAV{:});
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % equations of clock model
    dydt(1) = -1 * dm_per * Per + vmax_per * ( ...
        1 + fold_per * power(CLOCK_BMAL1 / Ka_per_cb, hill_per_cb)) / ( ...
            1 + power(CLOCK_BMAL1 / Ka_per_cb, hill_per_cb) * ...
            (1 + power(PER_CRY / Ki_per_pc, hill_per_pc)));

    dydt(2) = -1 * dm_cry * Cry + 1 / ( ...
        1 + power(REV / Ki_cry_rev, hill_cry_rev)) * ( ...
            vmax_cry * ...
            (1 + fold_cry * power(CLOCK_BMAL1 / Ka_cry_cb, hill_cry_cb))) / ( ...
                1 + power(CLOCK_BMAL1 / Ka_cry_cb, hill_cry_cb) * ...
                (1 + power(PER_CRY / Ki_cry_pc, hill_cry_pc)));

    dydt(3) = -1 * dm_rev * Rev + vmax_rev * ( ...
        1 + fold_rev * power(CLOCK_BMAL1 / Ka_rev_cb, hill_rev_cb)) / ( ...
            1 + power(CLOCK_BMAL1 / Ka_rev_cb, hill_rev_cb) * ...
            (1 + power(PER_CRY / Ki_rev_pc, hill_rev_pc)));

    dydt(4) = -1 * dm_ror * Ror + vmax_ror * ( ...
        1 + fold_ror * power(CLOCK_BMAL1 / Ka_ror_cb, hill_ror_cb)) / ( ...
            1 + power(CLOCK_BMAL1 / Ka_ror_cb, hill_ror_cb) * ...
            (1 + power(PER_CRY / Ki_ror_pc, hill_ror_pc)));

    dydt(5) = -1 * dm_bmal * Bmal1 + vmax_bmal * ( ...
        1 + fold_bmal * power(ROR / Ka_bmal_ror, hill_bmal_ror)) / ( ...
            1 + power(REV / Ki_bmal_rev, hill_bmal_rev) + ...
            power(ROR / Ka_bmal_ror, hill_bmal_ror));

    dydt(6) = -1 * dp_per * PER + kp_per * Per - ( ...
        kass_pc * PER * CRY - kdiss_pc * PER_CRY);

    dydt(7) = -1 * dp_cry * CRY + kp_cry * Cry - ...
        (kass_pc * PER * CRY - kdiss_pc * PER_CRY);

    dydt(8) = -1 * dp_rev * REV + kp_rev * Rev;

    dydt(9) = -1 * dp_ror * ROR + kp_ror * Ror;

    dydt(10) = -1 * dp_bmal * BMAL1 + kp_bmal * Bmal1 - ...
        (kass_cb * BMAL1 - kdiss_cb * CLOCK_BMAL1);

    dydt(11) = kass_pc * PER * CRY - kdiss_pc * PER_CRY - d_pc * PER_CRY;

    dydt(12) = kass_cb * BMAL1 - kdiss_cb * CLOCK_BMAL1 - d_cb * CLOCK_BMAL1;
    
    % equations of clock-controlled IAV model
    if (H < 0) || (t < t_IAV)
        H = 0;
        dydt(13) = 0;
    else
        dydt(13) = -1 * beta * H * V;
    end

    if (If < 0) || (t < t_IAV)
        If = 0;
        dydt(14) = 0;
    else
        dydt(14) = beta * H * V - alpha_MI * M * If - ...
            alpha_KI * K * If - alpha_TI * T_E * If / (K_TI + If) - d_I * If;
    end

    if (V < 0) || (t < t_IAV)
        V = 0;
        dydt(15) = 0;
    else
        dydt(15) = gamma * If - alpha_MV * M * V;
    end

    if (M < 0) || (t < t_IAV)
        M = 0;
        dydt(16) = 0;
    else
        dydt(16) = FracNoInf(K_B, K_B + BMAL1) * c_IM * ...
            FracNoInf(RealRootPromise(If, n_IM), (RealRootPromise(K_IM, n_IM) + RealRootPromise(If, n_IM))) ...
            + FracNoInf(K_R, K_R + REV) * c_M * FracNoInf(M, (K_M + M)) - d_M * M;
    end

    if (K < 0) || (t < t_IAV)
        K = 0;
        dydt(17) = 0;
    else
        dydt(17) = c_IK * FracNoInf(RealRootPromise(If, n_IK), (RealRootPromise(K_IK, n_IK) + RealRootPromise(If, n_IK))) + ...
            c_MK * FracNoInf(RealRootPromise(M, n_MK), (RealRootPromise(K_MK, n_MK) + RealRootPromise(M, n_MK))) - d_K * K;
    end

    if (T < 0) || (t < t_IAV)
        T = 0;
        dydt(18) = 0;
    else
        dydt(18) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
            FracNoInf(RealRootPromise(M, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(M, n_MT))) * T;
    end

    if (T_E < 0) || (t < t_IAV)
        T_E = 0;
        dydt(19) = 0;
    else
        dydt(19) = c_MT * FracNoInf(RealRootPromise(M, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(M, n_MT))) * T - d_T * T_E;
    end
    
end