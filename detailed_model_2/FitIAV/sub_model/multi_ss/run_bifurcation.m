% load base parameters
par_base = dlmread('multi_ss_par_256.txt', '', [par_base_num-1 0 par_base_num-1 57]);

par_consider = (Vnum - 1) * 1e-2;
par_consider_idx = 21;

tic;

[multi_ss_bool, multi_ss] = if_multi_ss(par_base, par_consider_idx, par_consider);

file_name = strcat('bifurcation_par_num_', int2str(par_base_num), '_Vnum_', int2str(Vnum), '.txt');
file_multi_ss = fopen(file_name, 'w');

fprintf(file_multi_ss, '%d ', par_base_num);
fprintf(file_multi_ss, '%f ', par_consider);
fprintf(file_multi_ss, '%f ', multi_ss_bool);
fprintf(file_multi_ss, '%f ', reshape(multi_ss, 1, []));
fprintf(file_multi_ss, '\n');

fclose(file_multi_ss);
toc;

disp(['run time:', num2str(toc)]);