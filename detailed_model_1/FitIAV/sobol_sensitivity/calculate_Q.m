function [Q_V, Q_M, Q_K] = calculate_Q( par_IAV_ini, par_consider_ind, par_consider, y_ini )
    
    y_consider = solve_odes(par_IAV_ini, par_consider_ind, par_consider);
    
    Q_V = 0;
    diff_V = (y_ini(:, 4) - y_consider(:, 4)) .^ 2;
    norm_val_V = max(y_ini(:, 4)) .^ 2;
    if norm_val_V ~= 0
        Q_V = trapz(diff_V) ./ norm_val_V;
    end
    
    Q_M = 0;
    diff_M = (y_ini(:, 5) - y_consider(:, 5)) .^ 2;
    norm_val_M = max(y_ini(:, 5)) .^ 2;
    if norm_val_M ~= 0
        Q_M = trapz(diff_M) ./ norm_val_M;
    end
    
    Q_K = 0;
    diff_K = (y_ini(:, 11) - y_consider(:, 11)) .^ 2;
    norm_val_K = max(y_ini(:, 11)) .^ 2;
    if norm_val_K ~= 0
        Q_K = trapz(diff_K) ./ norm_val_K;
    end
    
end

