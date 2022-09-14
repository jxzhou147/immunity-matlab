% load base parameters
par_base = importdata('par_M_init.txt');
par_base = par_base.data;
% translate to par in the odes
log_par_ind = (1:6);
for i = log_par_ind
    par_base(i) = 10 .^ par_base(i);
end

% load hls parameters
par_offset = arg1;
par_set_num = 10;
par_consider_idx = (1:7);
par_consider = dlmread('lhs_M.txt', '', [par_offset 0 par_set_num-1 6]);

multi_ss = false(par_set_num, 1);

p = parpool(20);
tic;
for i = 1:20
    multi_ss(i) = if_multi_ss(par_base, par_consider_idx, par_consider(i, :));
end
delete(p);
toc;

disp(['run time:', num2str(toc)]);

% write multi_ss to file
file_name = strcat('multi_ss_', int2str(par_offset), '.txt');
file_multi_ss = fopen(file_name, 'w');
for i = 1:par_set_num
    fprintf(file_multi_ss, '%d ', par_offset + i);
    fprintf(file_multi_ss, '%f\n', multi_ss(i));
end