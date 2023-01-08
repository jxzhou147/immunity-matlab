function [ infla ] = inflammation( Mono, N, par_infla )
    %fit inflammation and weight change from I, M, Mono, N
    
    sigma_Mono = par_infla(1);
    sigma_N = par_infla(2);
    
    infla = sigma_Mono * Mono + sigma_N * N;
    
end

