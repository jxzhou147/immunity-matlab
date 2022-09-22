% load base parameters
par_base = importdata('par_base.txt');
par_base = par_base.data;
% translate to par in the odes
log_par_ind = [1:39 43:54];
for i = log_par_ind
    par_base(i) = 10 .^ par_base(i);
end

% load hls parameters
par_consider_idx = (1:58);
par_consider = dlmread('multi_ss_par.txt', '', [par_num-1 0 par_num-1 57]);

tic;

[multi_ss_bool, multi_ss] = if_multi_ss(par_base, par_consider_idx, par_consider);

file_name = strcat('multi_ss_', int2str(par_num), '.txt');
file_multi_ss = fopen(file_name, 'w');
fprintf(file_multi_ss, '%d ', par_num);
fprintf(file_multi_ss, '%f ', multi_ss_bool);
fprintf(file_multi_ss, '%f ', reshape(multi_ss, 1, []));
fprintf(file_multi_ss, '\n');
fclose(file_multi_ss);
toc;

disp(['run time:', num2str(toc)]);
