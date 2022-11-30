% LHS for parameters
clear;clc;
% load parameters
par_init = importdata('par_base_IAV_clock_values.txt');
par_init = par_init.data;

ratio = 0.1;
N = 100;  % parameter set number

par_bound = zeros(2, length(par_init));
par_bound(1, :) = par_init .* (1 - ratio);
par_bound(2, :) = par_init .* (1 + ratio);

lhs_rad = lhsdesign(N, length(par_bound));
par_lhs = par_bound(1, :) + lhs_rad .* (par_bound(2, :) - par_bound(1, :));

% write lhs parameters to file
file_par = fopen('lhs_par.txt', 'w');
for i = 1:N
    fprintf(file_par, '%.10f ', par_lhs(i, :));
    fprintf(file_par, '\n');
end
fclose(file_par);
