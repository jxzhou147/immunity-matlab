function [ infla ] = inflammation( M, Mono, N, par_infla )
    %fit inflammation and weight change from I, M, Mono, N
    
    sigma_M = par_infla(1);
    sigma_Mono = par_infla(2);
    sigma_N = par_infla(3);
    
    infla = sigma_Mono * Mono + sigma_N * (N - 15);
%     infla = sigma_M * M + sigma_Mono * Mono + sigma_N * N;
    
end

