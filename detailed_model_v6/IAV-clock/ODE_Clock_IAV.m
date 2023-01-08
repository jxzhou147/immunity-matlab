function dydt = ODE_Clock_IAV(t, y, par_clock, par_IAV, t_IAV)
    
    dydt = zeros(25, 1);    % 1:12 for clock ODEs and 13:25 for IAV ODEs
    
    % variables of clock_IAV model
    y_cell = num2cell(y);
    [Per, Cry, Rev, Ror, Bmal1, PER, CRY, REV, ROR, ...
        BMAL1, PER_CRY, CLOCK_BMAL1, ...
        H, I1, I2, V, Mono, N, IL1b, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
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
    [zeta, beta, k_I, gamma, eta, r, a_MI, a_NI, a_KI, a_TI, ...
        c_IL1b_Mono, c_IL1b_N, c_CCL2_Mono, c_CXCL5_N, c_M_IL1b, c_I_IL1b, ...
        c_M_IL10, c_M_CCL2, c_I_CCL2, c_N_CXCL5, c_I_CXCL5, c_I_K, c_M_T, ...
        K_V, K_IL1b_Mono, K_IL1b_N, K_CCL2_Mono, K_CXCL5_N, K_IL10_IL1b, K_IL10_CCL2, K_IL10_CXCL5, ...
        K_I_K, K_T, K_M_T, n_CCL2_Mono, n_CXCL5_N, n_I_K, n_M_T, ...
        d_I, d_V, d_Vs, d_Mono, d_N, d_IL1b, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, b_H, b_K, ...
        K_BMAL1_IL1b, K_BMAL1_CCL2, K_REV_V, K_REV_IL10, ...
        K_REV_CCL2, K_REV_CXCL5, c_IL1b_REV, K_IL1b_REV, n_IL1b_REV] ...
        = deal(par_IAV{:});
    
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

    dydt(8) = -1 * dp_rev * (1 + c_IL1b_REV * ...
        FracNoInf(RealRootPromise(IL1b, n_IL1b_REV), (RealRootPromise(K_IL1b_REV, n_IL1b_REV) + RealRootPromise(IL1b, n_IL1b_REV)))) * ...
        REV + kp_rev * Rev;
%     dydt(8) = -1 * dp_rev * REV + kp_rev * Rev;

    dydt(9) = -1 * dp_ror * ROR + kp_ror * Ror;

    dydt(10) = -1 * dp_bmal * BMAL1 + kp_bmal * Bmal1 - ...
        (kass_cb * BMAL1 - kdiss_cb * CLOCK_BMAL1);

    dydt(11) = kass_pc * PER * CRY - kdiss_pc * PER_CRY - d_pc * PER_CRY;

    dydt(12) = kass_cb * BMAL1 - kdiss_cb * CLOCK_BMAL1 - d_cb * CLOCK_BMAL1;
    
    % equations of clock-controlled IAV model
    if (t < t_IAV)
        V = 0;
    end
    dydt(13) = zeta * H * (1 - H / b_H) - beta * H * V;

    dydt(14) = beta * H * V - k_I * I1;
    
    dydt(15) = k_I * I1 - a_MI * Mono * I2 - a_NI * N * I2 - a_KI  * K * I2 - a_TI * T_E * I2 - d_I * I2;

%     dydt(16) = gamma * I2 - d_V * V - d_Vs * V / (K_V + V);
    dydt(16) = FracNoInf(K_REV_V, K_REV_V + REV) * gamma * I2 - T_E / (T_E + 150) *  d_V * V - d_Vs * V / (K_V + V);
  
    dydt(17) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * ...
        FracNoInf(RealRootPromise(CCL2, n_CCL2_Mono), RealRootPromise(K_CCL2_Mono, n_CCL2_Mono) + RealRootPromise(CCL2, n_CCL2_Mono)) - d_Mono * Mono;

    dydt(18) = c_CXCL5_N * (1 + c_IL1b_N * FracNoInf(IL1b, (K_IL1b_N + IL1b))) * ...
        FracNoInf(RealRootPromise(CXCL5, n_CXCL5_N), RealRootPromise(K_CXCL5_N, n_CXCL5_N) + RealRootPromise(CXCL5, n_CXCL5_N)) - d_N * N;

    dydt(19) = FracNoInf(K_BMAL1_IL1b, K_BMAL1_IL1b + BMAL1) * FracNoInf(K_IL10_IL1b, K_IL10_IL1b + IL10) * c_M_IL1b * Mono + ...
        c_I_IL1b * I2 - d_IL1b * IL1b;

    dydt(20) = FracNoInf(K_REV_IL10, K_REV_IL10 + REV) * c_M_IL10 * Mono - d_IL10 * IL10;

    dydt(21) = FracNoInf(K_BMAL1_CCL2, K_BMAL1_CCL2 + BMAL1) * FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * Mono + ...
        FracNoInf(K_REV_CCL2, K_REV_CCL2 + REV) * c_I_CCL2 * I2 - d_CCL2 * CCL2;

    dydt(22) = FracNoInf(K_REV_CXCL5, K_REV_CXCL5 + REV) * FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_N_CXCL5 * N + ...
        c_I_CXCL5 * I2 - d_CXCL5 * CXCL5;

    dydt(23) = c_I_K * FracNoInf(RealRootPromise(I2, n_I_K), RealRootPromise(K_I_K, n_I_K) + RealRootPromise(I2, n_I_K)) - d_K * (K - b_K);

    dydt(24) = eta * T * (1 - FracNoInf(T, K_T)) - c_M_T * ...
        FracNoInf(RealRootPromise(Mono, n_M_T), (RealRootPromise(K_M_T, n_M_T) + RealRootPromise(Mono, n_M_T))) * T;

    dydt(25) = r * c_M_T * FracNoInf(RealRootPromise(Mono, n_M_T), (RealRootPromise(K_M_T, n_M_T) + RealRootPromise(Mono, n_M_T))) * T - d_T * T_E;

end