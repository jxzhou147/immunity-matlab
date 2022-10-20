% LHS for initial values
clear;clc
% load initial values bounds
init_bound = importdata('init_bound.txt');
init_bound = init_bound.data;

N = 1000;

lhs_rad = lhsdesign(N, length(init_bound));
init_lhs = lhs_rad .* init_bound(:, 1)' + lhs_rad .* (init_bound(:, 2)' - init_bound(:, 1)');

% write lhs parameters to file
file_init = fopen('lhs_init.txt', 'w');
for i = 1:N
    fprintf(file_init, '%f ', init_lhs(i, :));
    fprintf(file_init, '\n');
end
fclose(file_init);
