There is a big bug in all lhs_par.m !!!

par_lhs = lhs_rad .* par_bound(:, 1)' + lhs_rad .* (par_bound(:, 2)' - par_bound(:, 1)');

should be changed with

par_lhs = par_bound(:, 1)' + lhs_rad .* (par_bound(:, 2)' - par_bound(:, 1)');