function dydt = ODE_Clock_IAV(t, y, par_clock, par_IAV, t_IAV)
    
    dydt = zeros(25, 1);
    
    % variables of clock_IAV model
    y_cell = num2cell(y);
    [Per, Cry, Rev, Ror, Bmal1, PER, CRY, REV, ROR, ...
        BMAL1, PER_CRY, CLOCK_BMAL1, ...
        H, If, D, V, M, N, IL6, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
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
    [beta, gamma, eta, ...
        a_NH, a_MI, a_NI, a_KI, a_TI, a_NV, ...
        c_IL6_M, c_CCL2_M, c_CXCL5_N, c_M_IL6, c_K_IL6, c_N_IL6, c_I_IL6, ...
        c_IM, c_DM, c_M_IL10, c_M_CCL2, c_I_CCL2, c_M_CXCL5, c_I_CXCL5, c_IK, c_CCL2_K, c_MT, ...
        K_IL6_M, K_IM, K_DM, K_CCL2_M, K_CXCL5_N, K_IL10_IL6, K_IL10_CCL2, ...
        K_IL10_CXCL5, K_IK, K_CCL2_K, K_T, K_MT, ...
        n_IM, n_DM, n_CCL2_M, n_CXCL5_N, n_IK, n_CCL2_K, n_MT, ...
        d_H, d_I, d_V, d_M, d_N, d_IL6, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, ...
        b_H, b_M, b_N, b_K, ...
        K_B_M, K_BMAL1_IL6, c_P_K, K_P_K, K_REV_V, K_REV_IL6, K_REV_IL10, ...
        K_REV_CCL2, K_REV_CXCL5, c_ROR_IL6, K_ROR_IL6, c_IL6_REV, K_IL6_REV, n_IL6_REV] = deal(par_IAV{:});
    
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

    dydt(8) = -1 * dp_rev * (1 + c_IL6_REV * ...
        FracNoInf(RealRootPromise(IL6, n_IL6_REV), (RealRootPromise(K_IL6_REV, n_IL6_REV) + RealRootPromise(IL6, n_IL6_REV)))) * ...
        REV + kp_rev * Rev;

    dydt(9) = -1 * dp_ror * ROR + kp_ror * Ror;

    dydt(10) = -1 * dp_bmal * BMAL1 + kp_bmal * Bmal1 - ...
        (kass_cb * BMAL1 - kdiss_cb * CLOCK_BMAL1);

    dydt(11) = kass_pc * PER * CRY - kdiss_pc * PER_CRY - d_pc * PER_CRY;

    dydt(12) = kass_cb * BMAL1 - kdiss_cb * CLOCK_BMAL1 - d_cb * CLOCK_BMAL1;
    

    % equations of clock-controlled IAV model
    if (t <t_IAV)
        for i = 13:25
            dydt(i) = 0;
        end
    else
        if (H < 0)
            H = 0;
            dydt(13) = 0;
        else
            dydt(13) = d_H * b_H - beta * H * V - a_NH * (N - b_N) * H - d_H * H;
        end
    

        if (If < 0)
            If = 0;
            dydt(14) = 0;
        else
            dydt(14) = beta * H * V - FracNoInf(K_B_M, (K_B_M + BMAL1)) * (1 + c_IL6_M * FracNoInf(IL6, (K_IL6_M + IL6))) * a_MI * M * If - ...
                a_NI * N * If - (1 + FracNoInf(c_P_K * PER, (K_P_K + PER))) * a_KI * K * If - a_TI * T_E * If - d_I * If;
        end

        if (D < 0)
            D = 0;
            dydt(15) = 0;
        else
            dydt(15) = FracNoInf(K_B_M, (K_B_M + BMAL1)) * (1 + c_IL6_M * FracNoInf(IL6, (K_IL6_M + IL6))) * a_MI * M * If + ...
                a_NI * N * If + (1 + FracNoInf(c_P_K * PER, (K_P_K + PER))) * a_KI * K * If + a_TI * T_E * If + d_I * If + ...
                a_NH * (N - b_N) * H + d_H * H;
        end

        if (V < 0)
            V = 0;
            dydt(16) = 0;
        else
            dydt(16) = FracNoInf(K_REV_V, (K_REV_V + REV)) * gamma * If - a_NV * N * V - d_V * V;
        end

        if (M < 0)
            M = 0;
            dydt(17) = 0;
        else
            dydt(17) = b_M * d_M + c_IM * FracNoInf(RealRootPromise(If, n_IM), (RealRootPromise(K_IM, n_IM) + RealRootPromise(If, n_IM))) + ...
                c_DM * FracNoInf(RealRootPromise(D, n_DM), (RealRootPromise(K_DM, n_DM) + RealRootPromise(D, n_DM))) + ...
                c_CCL2_M * FracNoInf(RealRootPromise(CCL2, n_CCL2_M), (RealRootPromise(K_CCL2_M, n_CCL2_M) + RealRootPromise(CCL2, n_CCL2_M))) - ...
                d_M * M;
        end

        if (N < 0)
            N = 0;
            dydt(18) = 0;
        else
            dydt(18) = b_N * d_N + c_CXCL5_N * FracNoInf(RealRootPromise(CXCL5, n_CXCL5_N), (RealRootPromise(K_CXCL5_N, n_CXCL5_N) + RealRootPromise(CXCL5, n_CXCL5_N))) - ...
                d_N * N;
        end

        if (IL6 < 0)
            IL6 = 0;
            dydt(19) = 0;
        else
            dydt(19) =  FracNoInf(K_BMAL1_IL6, K_BMAL1_IL6 + BMAL1) * FracNoInf(K_REV_IL6, K_REV_IL6 +REV) * ...
                (1 + c_ROR_IL6 * FracNoInf(ROR, K_ROR_IL6 + ROR)) * FracNoInf(K_IL10_IL6, K_IL10_IL6 + IL10) * c_M_IL6 * M + ...
                c_K_IL6 * K + c_N_IL6 * N + c_I_IL6 * If - d_IL6 * IL6;
        end

        if (IL10 < 0)
            IL10 = 0;
            dydt(20) = 0;
        else
            dydt(20) = FracNoInf(K_REV_IL10, K_REV_IL10 + REV) * c_M_IL10 * M - d_IL10 * IL10;
        end

        if (CCL2 < 0)
            CCL2 = 0;
            dydt(21) = 0;
        else
            dydt(21) = FracNoInf(K_REV_CCL2, K_REV_CCL2 + REV) * FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * M + ...
                c_I_CCL2 * If - d_CCL2 * CCL2;
        end

        if (CXCL5 < 0)
            CXCL5 = 0;
            dydt(22) = 0;
        else
            dydt(22) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_M_CXCL5 * M + ...
                FracNoInf(K_REV_CXCL5, K_REV_CXCL5 + REV) * c_I_CXCL5 * If - d_CXCL5 * CXCL5;
        end

        if (K < 0)
            K = 0;
            dydt(23) = 0;
        else
            dydt(23) = b_K * d_K + c_IK * FracNoInf(RealRootPromise(If, n_IK), (RealRootPromise(K_IK, n_IK) + RealRootPromise(If, n_IK))) + ...
                c_CCL2_K * FracNoInf(RealRootPromise(CCL2, n_CCL2_K), RealRootPromise(K_CCL2_K, n_CCL2_K) + RealRootPromise(CCL2, n_CCL2_K)) ...
                - d_K * K;
        end

        if (T < 0)
            T = 0;
            dydt(24) = 0;
        else
            dydt(24) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
                FracNoInf(RealRootPromise(M, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(M, n_MT))) * T;
        end

        if (T_E < 0)
            T_E = 0;
            dydt(25) = 0;
        else
            dydt(25) = c_MT * FracNoInf(RealRootPromise(M, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(M, n_MT))) * T - d_T * T_E;
        end
    end
    
end