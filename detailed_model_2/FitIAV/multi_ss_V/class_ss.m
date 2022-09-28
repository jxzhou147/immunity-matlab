multi_ss = importdata('merged_multi_ss.txt');

file_id = fopen('class_merged_multi_ss.txt', 'w');
for i = 3:2:30
    diff_idx = multi_ss(multi_ss(:, i) ~= multi_ss(:, i+1));
    fprintf(file_id, '%d ', diff_idx);
    fprintf(file_id, '\n');
end