% LHS for initial values
var_name = ['Mono', 'CCL2'];
var_num = 2;

var_bound = zeros(2, var_num);
var_bound(:, 1) = [0; 1e3];
var_bound(:, 2) = [0; 1e4];

N = 1000;

lhs_rad = lhsdesign(N, var_num);
var_lhs = lhs_rad .* var_bound(1, :) + lhs_rad .* (var_bound(2, :) - var_bound(1, :));

% write lhs parameters to file
file_var_M = fopen('lhs_init.txt', 'w');
for i = 1:N
    fprintf(file_var_M, '%f ', var_lhs(i, :));
    fprintf(file_var_M, '\n');
end
fclose(file_var_M);
