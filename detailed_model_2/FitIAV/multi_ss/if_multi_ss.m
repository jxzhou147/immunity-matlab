function multi_ss = if_multi_ss(par_base, par_consider_idx, par_consider)
    
    multi_ss = false;
    
    % load initial values
    y0s = importdata('lhs_init.txt');
    [m, n] = size(y0s);
    
    ys = zeros(m, n);
    
    parfor i = 1:m
        ys(i, :) = solve_odes(par_base, par_consider_idx, par_consider, y0s(i, :));
    end
    
    for i = 1:n
        if ((max(ys(:, i)) - min(ys(:, i))) / max(ys(:, i)) > 1e-1) && (max(ys(:, i)) - min(ys(:, i)) > 1)
            multi_ss = true;
        end
    end
    
end