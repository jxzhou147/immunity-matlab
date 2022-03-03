function dydt = ODE_Clock(t, y, par)
    
    dydt = zeros(12);
    
    % variables of clock model
    y_cell = num2cell(y);
    [Per, Cry, Rev, Ror, Bmal1, PER, CRY, REV, ROR, ...
        BMAL1, PER_CRY, CLOCK_BMAL1] = deal(y_cell{:});
    
    % parameters
    % mRNA and protein degradation rate constants (in h^-1)
    par1 = num2cell(par(1:12));
    [dm_per, dm_cry, dm_rev, dm_ror, dm_bmal, dp_per, ...
        dp_cry, dp_rev, dp_ror, dp_bmal, d_pc, d_cb] = deal(par1{:});
    
    % maximal transcription rates (in nmol l^-1 h^-1)
    par2 = num2cell(par(13:17));
    [vmax_per, vmax_cry, vmax_rev, vmax_ror, vmax_bmal] = deal(par2{:});
    
    % activation ratios (demensionless)
    par3 = num2cell(par(18:22));
    [fold_per, fold_cry, fold_rev, fold_ror, fold_bmal] = deal(par3{:});
    
    % regulation thresholds (in nmol/l)
    par4 = num2cell(par(23:33));
    [Ka_per_cb, Ki_per_pc, Ka_cry_cb, Ki_cry_pc, Ki_cry_rev, Ka_rev_cb, ...
        Ki_rev_pc, Ka_ror_cb, Ki_ror_pc, Ka_bmal_ror, Ki_bmal_rev] = deal(par4{:});
    
    % hill coefficients (dimensionless)
    par5 = num2cell(par(34:44));
    [hill_per_cb, hill_per_pc, hill_cry_cb, hill_cry_pc, ...
        hill_cry_rev, hill_rev_cb, hill_rev_pc, hill_ror_cb, ...
        hill_ror_pc, hill_bmal_ror, hill_bmal_rev] = deal(par5{:});
    
    % translation rates (in molecules per hour per mRNA)
    par6 = num2cell(par(45:49));
    [kp_per, kp_cry, kp_rev, kp_ror, kp_bmal] = deal(par6{:});
    
    % complexation kinetic rates
    par7 = num2cell(par(50:53));
    [kass_cb, kass_pc, kdiss_cb, kdiss_pc] = deal(par7{:});