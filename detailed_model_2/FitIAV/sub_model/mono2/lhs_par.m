% LHS for parameters
par_name = ['c_CCL2_Mono', 'K_CCL2_Mono', 'n_CCL2_Mono', 'd_Mono', 'c_M_CCL2', 'c_I_CCL2', 'd_CCL2', 'I'];

% load initial parameters
par_init = importdata('par_M_init.txt');
par_init = par_init.data;

% translate to par in the odes
log_par_ind = (1:7);
for i = log_par_ind
    par_init(i) = 10 .^ par_init(i);
end

par_consider = par_init(1:8);
par_consider_num = length(par_consider);

% ratio to define parameter bounds and sample number
ratio = 10;
N = 1000;

par_bound = zeros(2, par_consider_num);
par_bound(1, :) = par_consider ./ ratio;
par_bound(2, :) = par_consider .* ratio;
par_bound(1, 8) = 0;
par_bound(2, 8) = 5e3;

lhs_rad = lhsdesign(N, par_consider_num);
par_lhs = lhs_rad .* par_bound(1, :) + lhs_rad .* (par_bound(2, :) - par_bound(1, :));

% write lhs parameters to file
file_par_M = fopen('lhs_M.txt', 'w');
for i = 1:N
    fprintf(file_par_M, '%f ', par_lhs(i, :));
    fprintf(file_par_M, '\n');
end
fclose(file_par_M);
