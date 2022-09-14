function dydt = ODE_M(t, y, par_M)
    
    dydt = zeros(2, 1);
    
    Mono = y(1);
    CCL2 = y(2);
    
    par_M = num2cell(par_M);
    [c_CCL2_Mono, Mono_tot, d_Mono, c_M_CCL2, c_I_CCL2, d_CCL2, I] = deal(par_M{:});
    
    dydt(1) = c_CCL2_Mono * CCL2 * (Mono_tot - Mono) - d_Mono * Mono;
    
    dydt(2) = c_M_CCL2 * Mono + c_I_CCL2 * I - d_CCL2 * CCL2;
    
end