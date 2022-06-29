function dydt = ODE_Clock_IAV(t, y, par_clock, par_IAV, t_IAV)
    
    dydt = zeros(27, 1);
    
    % variables of clock_IAV model
    y_cell = num2cell(y);
    [Per, Cry, Rev, Ror, Bmal1, PER, CRY, REV, ROR, ...
        BMAL1, PER_CRY, CLOCK_BMAL1, ...
        H, If, D, V, M_0, M, Mono, N, IL1b, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
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
        a_NH, a_MH, a_MI, a_NI, a_KI, a_TI, a_MD, a_NV, ...
        c_IL1b_M, c_IL1b_Mono, c_IL1b_N, c_IL1b_K, c_IM, c_CCL2_Mono, c_N, c_CXCL5_N, c_M_IL1b, c_I_IL1b, c_D_IL1b, ...
        c_M_IL10, c_M_CCL2, c_I_CCL2, c_N_CXCL5, c_I_CXCL5, c_K, c_MK, c_IK, c_MT, ...
        K_IL1b_M, K_IL1b_Mono, K_IL1b_N, K_IL1b_K, K_MI, K_I_M, K_D_IL1b, K_IL10_IL1b, K_IL10_CCL2, K_IL10_CXCL5, K_T, K_MT, n_I_M, n_D_IL1b, n_MT, ...
        d_H, d_I, d_V, d_M0, d_M, d_Mono, d_N, d_IL1b, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, ...
        Mono_tot, b_H, b_M0, b_N, b_K, ...
        K_B_M, K_BMAL1_IL1b, c_P_K, K_P_K, K_REV_V, K_REV_IL1b, K_REV_IL10, ...
        K_REV_CCL2, K_REV_CXCL5, c_ROR_IL1b, K_ROR_IL1b, c_IL1b_REV, K_IL1b_REV, n_IL1b_REV] ...
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
    if (t >= t_IAV)
        dydt(13) = (a_NH * b_N + a_MH * b_M0 + d_H) * b_H - beta * H * V - ...
            a_NH * N * H - a_MH * (M + Mono) * H - d_H * H;
        
        dydt(14) = beta * H * V - FracNoInf(K_B_M, (K_B_M + BMAL1)) * a_MI * FracNoInf((M + Mono), (K_MI + M + Mono)) * If - ...
            a_NI * N * If - (1 + FracNoInf(c_P_K * PER, (K_P_K + PER))) * a_KI * K * If - a_TI * T_E * If - d_I * If;

        dydt(15) = a_NI * N * If + (1 + FracNoInf(c_P_K * PER, (K_P_K + PER))) * a_KI * K * If + a_TI * T_E * If + d_I * If + ...
            a_NH * N * H + a_MH * (M + Mono) * H - a_MD * (M + Mono) * D;

        dydt(16) = FracNoInf(K_REV_V, (K_REV_V + REV)) * gamma * If - a_NV * N * V - d_V * V;
        
        dydt(17) = b_M0 * d_M0 - (1 + c_IL1b_M * FracNoInf(IL1b, K_IL1b_M +IL1b)) * ...
            c_IM * FracNoInf(RealRootPromise(If, n_I_M), RealRootPromise(K_I_M, n_I_M) + RealRootPromise(If, n_I_M)) * M_0 - ...
            d_M0 * M_0;
        
        dydt(18) = (1 + c_IL1b_M * FracNoInf(IL1b, K_IL1b_M +IL1b)) * ...
            c_IM * FracNoInf(RealRootPromise(If, n_I_M), RealRootPromise(K_I_M, n_I_M) + RealRootPromise(If, n_I_M)) * M_0 - ...
            d_M * M;
        
        dydt(19) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * 1000 * FracNoInf(CCL2 .^ 3, CCL2 .^ 3 + 2000 .^ 3) * (Mono_tot - Mono) - d_Mono * Mono;

        dydt(20) = (c_N + c_CXCL5_N * (1 + c_IL1b_N * FracNoInf(IL1b, (K_IL1b_N + IL1b))) * CXCL5) * ((c_N + d_N) / c_N * b_N - N) - d_N * N;

%         dydt(21) = (FracNoInf(K_BMAL1_IL1b, K_BMAL1_IL1b + BMAL1) * FracNoInf(K_REV_IL1b, K_REV_IL1b +REV) * ...
%             (1 + c_ROR_IL1b * FracNoInf(ROR, K_ROR_IL1b + ROR)) * FracNoInf(K_IL10_IL1b, K_IL10_IL1b + IL10) * c_M_IL1b * (M + Mono) + FracNoInf(K_REV_IL1b, K_REV_IL1b +REV) * c_I_IL1b * If) * ...
%             (1 + c_D_IL1b * FracNoInf(RealRootPromise(D, n_D_IL1b), RealRootPromise(K_D_IL1b, n_D_IL1b) + RealRootPromise(D, n_D_IL1b))) - d_IL1b * IL1b;

        dydt(21) = (FracNoInf(K_BMAL1_IL1b, K_BMAL1_IL1b + BMAL1) * FracNoInf(K_REV_IL1b, K_REV_IL1b +REV) * ...
            (1 + c_ROR_IL1b * FracNoInf(ROR, K_ROR_IL1b + ROR)) * FracNoInf(K_IL10_IL1b, K_IL10_IL1b + IL10) * c_M_IL1b * (M + Mono) + FracNoInf(K_REV_IL1b, K_REV_IL1b +REV) * c_I_IL1b * (If + D)) - d_IL1b * IL1b;

        dydt(22) = FracNoInf(K_REV_IL10, K_REV_IL10 + REV) * c_M_IL10 * (M + Mono) - d_IL10 * IL10;

%         if ((t_IAV -115) < 1)
%             if (t > (t_IAV + 10))
%                 REV = 0;
%             end
% %         elseif ((t_IAV -127) < 1)
% %             if (t > (t_IAV + 30))
% %                 REV = 0;
% %             end
%         end
    %     c_M_CCL2 = (1 - exp(-0.1 * t)) * c_M_CCL2;
        
        dydt(23) = FracNoInf(K_REV_CCL2, K_REV_CCL2 + REV) * FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * (M + Mono) + ...
            c_I_CCL2 * If - d_CCL2 * CCL2;

        dydt(24) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_N_CXCL5 * N + ...
            FracNoInf(K_REV_CXCL5, K_REV_CXCL5 + REV) * c_I_CXCL5 * If - d_CXCL5 * CXCL5;

        dydt(25) = (1 + c_IL1b_K * FracNoInf(IL1b, (K_IL1b_K + IL1b))) * (c_K + c_MK * M + c_IK * If) * ((c_K + d_K) / c_K * b_K - K) - d_K * K;

    %    dydt(13) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
    %         FracNoInf(RealRootPromise((M + Mono), n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise((M + Mono), n_MT))) * T;
    % 
    %     dydt(14) = c_MT * FracNoInf(RealRootPromise((M + Mono), n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise((M + Mono), n_MT))) * T - d_T * T_E;
        dydt(26) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
            FracNoInf(RealRootPromise(K, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(K, n_MT))) * T;

        dydt(27) = c_MT * FracNoInf(RealRootPromise(K, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(K, n_MT))) * T - d_T * T_E;
    end
end