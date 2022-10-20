function dydt = ODE_IAV(t, y, par_IAV)

%     for i = 1:length(y)
%         if y(i) < 1
%             y(i) = 0;
%         end
%     end
    
    dydt = zeros(14, 1);    % y(16):weight change
    
    % variables of IAV model
    y_cell = num2cell(y);
    [H, If, V, M_0, M, Mono, N, IL1b, IL10, CCL2, CXCL5, K, T, T_E] = deal(y_cell{:});
    
    % parameters of IAV model
    par_IAV = num2cell(par_IAV);
    [zeta, beta, gamma, eta, ...
        a_MI, a_NI, a_KI, a_TI, a_NV, ...
        c_IL1b_M, c_IL1b_Mono, c_IL1b_N, c_IM, c_CCL2_Mono, c_N, c_CXCL5_N, c_M_IL1b, c_I_IL1b, ...
        c_M_IL10, c_M_CCL2, c_I_CCL2, c_N_CXCL5, c_I_CXCL5, c_K, c_MK, c_IK, c_MT, ...
        K_IL1b_M, K_IL1b_Mono, K_IL1b_N, K_MI, K_I_M, K_CCL2_Mono, K_IL10_IL1b, K_IL10_CCL2, K_IL10_CXCL5, K_T, K_MT, n_I_M, n_CCL2_Mono, n_MT, ...
        d_I, d_V, d_M0, d_M, d_Mono, d_N, d_IL1b, d_IL10, d_CCL2, d_CXCL5, d_K, d_T, ...
        b_H, b_M0, b_N, b_K] = deal(par_IAV{:});
    
    % parameters of inflammation and weigth change
%     par_infla = num2cell(par_infla);
%     [sigma_M, sigma_Mono, sigma_N, sigma_I, sigma_Infla, K_I, K_In, n_I, n_In, d_wc] = deal(par_infla{:});
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % equations of clock-controlled IAV model
%     dydt(1) = (a_NH * b_N + a_MH * b_M0 + d_H) * b_H - beta * H * V - ...
%         a_NH * N * H - a_MH * (M + Mono) * H - d_H * H;
%     dydt(1) = d_H * b_H - beta * H * V - a_NH * (N - b_N) * H - d_H * H;
%     dydt(1) = - beta * H * V;
%     if (y(1) < 0)
%         y(1) = 0;
%         dydt(1) = 0;
%     else
    
    dydt(1) = zeta * H * (1 - H / b_H ) - beta * H * V;
%     end
    
    dydt(2) = beta * H * V - a_MI * FracNoInf((M + Mono), (K_MI + M + Mono)) * If - ...
        a_NI * N * If - a_KI * K * If - a_TI * T_E * If - d_I * If;
%     dydt(2) = beta * H * V - a_MI * FracNoInf((M + Mono), (K_MI + M + Mono)) * If - ...
%         a_NI * N * If - a_KI * (1 - FracNoInf(RealRootPromise(M + Mono, 5), 50 ^ 5 + RealRootPromise(M + Mono, 5))) * K  * If - a_TI * (1 - FracNoInf(RealRootPromise(M + Mono, 5), 50 ^ 5 + RealRootPromise(M + Mono, 5))) * T_E * If * 200 / (200 + If) - d_I * If;
%     dydt(2) = beta * H * V - a_MI * FracNoInf((M + Mono), (K_MI + M + Mono)) * If - ...
%         a_NI * N * If - a_KI * FracNoInf(70 ^ 5, 70 ^ 5 + RealRootPromise(M + Mono, 5)) * K  * If - a_TI * FracNoInf(70 ^ 5, 70 ^ 5 + RealRootPromise(M + Mono, 5)) * T_E * If * 200 / (200 + If) - d_I * If;
   
%     dydt(3) = a_NI * N * If + a_KI * K * If + a_TI * T_E * If + d_I * If + ...
%         a_NH * N * H + a_MH * (M + Mono) * H - a_MD * (M + Mono) * D;

%     if (y(3) <1e-3)
%         dydt(3) = 0;
%     else
%     if V < 1
%         V = 0;
%     else
%     dydt(3) = gamma * If - a_NV * N * V - d_V * V;
%     end
%     end
%     y(3) = V;
%     dydt(3) = 0;
    
    % the last term is used to avoid regeneration of virus which is tiny enough
%     dydt(3) = gamma * If - a_NV * N * V - d_V * V - 8e-2 * V / (1e-4 + V);
    dydt(3) = gamma * If - a_NV * N * V - d_V * V - 8e-2 * V / (1e-2 + V);

    dydt(4) = b_M0 * d_M0 - (1 + c_IL1b_M * FracNoInf(IL1b, K_IL1b_M +IL1b)) * ...
        c_IM * FracNoInf(RealRootPromise(If, n_I_M), RealRootPromise(K_I_M, n_I_M) + RealRootPromise(If, n_I_M)) * M_0 - ...
        d_M0 * M_0;

    dydt(5) = (1 + c_IL1b_M * FracNoInf(IL1b, K_IL1b_M +IL1b)) * ...
        c_IM * FracNoInf(RealRootPromise(If, n_I_M), RealRootPromise(K_I_M, n_I_M) + RealRootPromise(If, n_I_M)) * M_0 - ...
        d_M * M;

%     dydt(7) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * CCL2 * (Mono_tot - Mono) - d_Mono * Mono;
%     dydt(7) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * 1000 * FracNoInf(CCL2, CCL2 + 500) * (Mono_tot - Mono) - d_Mono * Mono;
    dydt(6) = c_CCL2_Mono * (1 + c_IL1b_Mono * FracNoInf(IL1b, (K_IL1b_Mono + IL1b))) * ...
        FracNoInf(RealRootPromise(CCL2, n_CCL2_Mono), RealRootPromise(K_CCL2_Mono, n_CCL2_Mono) + RealRootPromise(CCL2, n_CCL2_Mono)) - d_Mono * Mono;
    
    dydt(7) = (c_N + c_CXCL5_N * (1 + c_IL1b_N * FracNoInf(IL1b, (K_IL1b_N + IL1b))) * CXCL5) * ((c_N + d_N) / c_N * b_N - N) - d_N * N;

%     dydt(9) = (FracNoInf(K_IL10_IL1b, K_IL10_IL1b + IL10) * c_M_IL1b * (M + Mono) + c_I_IL1b * If) * ...
%         (1 + c_D_IL1b * FracNoInf(RealRootPromise(D, n_D_IL1b), RealRootPromise(K_D_IL1b, n_D_IL1b) + RealRootPromise(D, n_D_IL1b))) - d_IL1b * IL1b;
    dydt(8) = (FracNoInf(K_IL10_IL1b, K_IL10_IL1b + IL10) * c_M_IL1b * (M + Mono) + c_I_IL1b * If) - d_IL1b * IL1b;

    dydt(9) = c_M_IL10 * (M + Mono) - d_IL10 * IL10;

%     c_M_CCL2 = (1 - exp(-0.1 * t)) * c_M_CCL2;
    dydt(10) = FracNoInf(K_IL10_CCL2, K_IL10_CCL2 + IL10) * c_M_CCL2 * (M + Mono) + ...
        c_I_CCL2 * If - d_CCL2 * CCL2;

%     dydt(11) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_N_CXCL5 * N + ...
%         c_I_CXCL5 * If - d_CXCL5 * CXCL5;
    dydt(11) = FracNoInf(K_IL10_CXCL5, K_IL10_CXCL5 + IL10) * c_N_CXCL5 * (N - b_N) + ...
        c_I_CXCL5 * If - d_CXCL5 * CXCL5;
    
    dydt(12) = (c_K + c_MK * M + c_IK * If) * ((c_K + d_K) / c_K * b_K - K) - d_K * K;
%     dydt(13) = (1 + c_IL1b_K * FracNoInf(IL1b, (K_IL1b_K + IL1b))) * (c_K + c_MK * M + c_IK * 20 * FracNoInf(If, 10 + If)) * ((c_K + d_K) / c_K * b_K - K) - d_K * K;

%    dydt(13) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
%         FracNoInf(RealRootPromise((M + Mono), n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise((M + Mono), n_MT))) * T;
% 
%     dydt(14) = c_MT * FracNoInf(RealRootPromise((M + Mono), n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise((M + Mono), n_MT))) * T - d_T * T_E;
    dydt(13) = eta * T * (1 - FracNoInf(T, K_T)) - c_MT * ...
        FracNoInf(RealRootPromise(Mono, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(Mono, n_MT))) * T;

    dydt(14) = c_MT * FracNoInf(RealRootPromise(Mono, n_MT), (RealRootPromise(K_MT, n_MT) + RealRootPromise(Mono, n_MT))) * T - d_T * T_E;
    
    %%%%%%%%%%%%%%%% dynamics of inflammation and weigth change
%     infla = sigma_M * M + sigma_Mono * Mono + sigma_N * N;
%     infla = sigma_Mono * Mono + sigma_N * (N - b_N);
%     dydt(16) = sigma_I * FracNoInf(RealRootPromise(If, n_I), RealRootPromise(K_I, n_I) + RealRootPromise(If, n_I)) + ...
%             sigma_Infla * FracNoInf(RealRootPromise(infla, n_In), RealRootPromise(K_In, n_In) + RealRootPromise(infla, n_In)) - d_wc * wc;
    

end