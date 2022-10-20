% load base parameters
par_base = importdata('par_M_init.txt');
par_base = par_base.data;
% translate to par in the odes
log_par_ind = (1:7);
for i = log_par_ind
    par_base(i) = 10 .^ par_base(i);
end

% load hls parameters
par_set_num = 100;
par_consider_idx = (1:8);
par_consider = dlmread('lhs_M.txt', '', [par_offset 0 par_set_num-1 7]);

multi_ss_bool = false(par_set_num, 1);
multi_ss = zeros(par_set_num, 2);

p = parpool(20);
tic;
for i = 1:par_set_num
    [multi_ss_bool(i), multi_ss(i, :)] = if_multi_ss(par_base, par_consider_idx, par_consider(i, :));
end
delete(p);
toc;

disp(['run time:', num2str(toc)]);

% write multi_ss to file
file_name = strcat('multi_ss_', int2str(par_offset), '.txt');
file_multi_ss = fopen(file_name, 'w');
for i = 1:par_set_num
    fprintf(file_multi_ss, '%d ', par_offset + i);
    fprintf(file_multi_ss, '%f ', multi_ss_bool(i));
    fprintf(file_multi_ss, '%f ', multi_ss(i, :));
    fprintf(file_multi_ss, '\n');
end