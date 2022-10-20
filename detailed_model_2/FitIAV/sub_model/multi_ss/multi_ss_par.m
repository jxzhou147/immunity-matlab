% find the multi_ss parameters
multi_ss_par_idx = [];

for i = 0:0
%     file_name = strcat('multi_ss_', int2str(i*100), '.txt');
    file_name = 'merged_multi_ss.txt';
    data_i = importdata(file_name);
    non_zero_idx = find(data_i(:, 2));
    multi_ss_par_idx_i = data_i(non_zero_idx, 1);
    multi_ss_par_idx = [multi_ss_par_idx; multi_ss_par_idx_i];
end

par_set = importdata('lhs_par.txt');
file_multi_ss_par = fopen('multi_ss_par.txt', 'w');
for i = 1:length(multi_ss_par_idx)
    fprintf(file_multi_ss_par, '%f ', par_set(multi_ss_par_idx(i), :));
    fprintf(file_multi_ss_par, '\n');
end
fclose(file_multi_ss_par);
