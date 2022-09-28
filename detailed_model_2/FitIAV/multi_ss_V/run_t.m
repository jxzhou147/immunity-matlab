% load parameter set
par_base = importdata('multi_ss_par_488.txt');
par_base = par_base(par_num, :);
par_consider_idx = (1:57);

% load initial values
y0s = importdata('lhs_init.txt');
y0s = y0s(1:1000, :);
[m, n] = size(y0s);

tspan = 0:1:10000;

y = zeros(length(tspan), n, m);

p = parpool(20);
tic;
parfor i = 1:m
    y_i = solve_odes(par_base, par_consider_idx, par_base, y0s(i, :), tspan);
    if size(y_i) == [length(tspan) n]
        y(:, :, i) = y_i;
    else
        y(:, :, i) = 404;
    end
end
toc;
disp(['run time:', num2str(toc)]);
delete(p);

% write curves to file
file_name = strcat('curve_parnum_', int2str(par_num), '.txt');
file_id = fopen(file_name, 'w');
for i = 1:m
    for j = 1:n
        fprintf(file_id, '%f ', y(:, j, i));
        fprintf(file_id, '\n');
    end
end
fclose(file_id);