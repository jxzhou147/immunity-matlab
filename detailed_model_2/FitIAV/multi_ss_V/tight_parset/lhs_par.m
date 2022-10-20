% LHS for parameters
clear;clc;
% load parameter bounds
par_bound = importdata('par_bound.txt');
par_bound = par_bound.data;

% translate to par in the odes
log_par_ind = [1:38 42:53];
for i = log_par_ind
    par_bound(i, :) = 10 .^ par_bound(i, :);
end

par_consider_idx = [14, 20, 21, 33, 40, 46, 50];
par_bound = par_bound(par_consider_idx, :);

N = 1000;  % parameter set number

lhs_rad = lhsdesign(N, length(par_bound));
par_lhs = lhs_rad .* par_bound(:, 1)' + lhs_rad .* (par_bound(:, 2)' - par_bound(:, 1)');

% write lhs parameters to file
file_par = fopen('lhs_par.txt', 'w');
for i = 1:N
    fprintf(file_par, '%.10f ', par_lhs(i, :));
    fprintf(file_par, '\n');
end
fclose(file_par);
