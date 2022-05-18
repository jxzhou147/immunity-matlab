clear;clc

% load initial IAV parameters
par_IAV_ini = importdata('par_IAV.txt');
par_IAV_ini = par_IAV_ini.data;

% indices of analysed parameters
par_consider_ind = [5, 7, 2, 16, 20, 21, 24];

% load generated parameters
par_consider = load('par_values.txt');

% solves odes for initial parameters
par_consider_ini = par_IAV_ini(par_consider_ind);
y_ini = solve_odes(par_IAV_ini, par_consider_ind, par_consider_ini);

% generate Q values for each parameters using parfor
tic;

par_num = length(par_consider);

tmax = 250;
tspan = 0:1:tmax;
y_Vs = zeros(par_num, tmax +1);
y_Ms = zeros(par_num, tmax +1);
y_Ks = zeros(par_num, tmax + 1);

p = parpool(20);

parfor i = 1:par_num
    y = solve_odes(par_IAV_ini, par_consider_ind, par_consider(i, :));
    
    y_Vs(i, :) = y(:, 4);
    y_Ms(i, :) = y(:, 5);
    y_Ks(i, :) = y(:, 11);
end

delete(p);

% write y to files
file_yV = fopen('y_V.txt', 'w');
for i = 1:par_num
    fprintf(file_yV, '%f ', y_Vs(i, :));
    fprintf(file_yV, '\n');
end
fclose(file_yV);

file_yM = fopen('y_M.txt', 'w');
for i = 1:par_num
    fprintf(file_yM, '%f ', y_Ms(i, :));
    fprintf(file_yM, '\n');
end
fclose(file_yM);

file_yK = fopen('y_K.txt', 'w');
for i = 1:par_num
    fprintf(file_yK, '%f ', y_Ks(i, :));
    fprintf(file_yK, '\n');
end
fclose(file_yK);

toc;
disp(['run time:', num2str(toc)]);
