function dydt = ODE_Mt(t, y, par_M, t_IAV)
    
    dydt = zeros(2, 1);
    
    Mono = y(1);
    CCL2 = y(2);
    
    par_M = num2cell(par_M);
    [c_CCL2_Mono, K_CCL2_Mono, n_CCL2_Mono, d_Mono, c_M_CCL2, c_I_CCL2, d_CCL2, I] = deal(par_M{:});
     
    c_I_CCL2 = 1 * (sin(2 * pi / 24 * t - pi / 2) + 1) + 1;
%     if (t >= t_IAV + 10)
%         c_I_CCL2 = 0.3;
%     end
    
    if (t >= t_IAV)
        dydt(1) = c_CCL2_Mono * CCL2 ^ n_CCL2_Mono / (K_CCL2_Mono ^ n_CCL2_Mono + CCL2 ^ n_CCL2_Mono) - d_Mono * Mono;

        dydt(2) = c_M_CCL2 * Mono + c_I_CCL2 * I - d_CCL2 * CCL2;
    end
    
end