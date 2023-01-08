function [ infla ] = inflammation( Mono, N, par_infla )
    %fit inflammation and weight change from I, M, Mono, N
    
    sigma_Mono = par_infla(2);
    sigma_N = par_infla(3);
    
    infla = sigma_Mono * Mono + sigma_N * N;
    
end

