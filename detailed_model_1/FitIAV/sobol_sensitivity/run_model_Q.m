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
QVs = zeros(par_num, 1);
QMs = zeros(par_num, 1);
QKs = zeros(par_num, 1);

p = parpool(20);

parfor i = 1:par_num
    [QVs(i), QMs(i), QKs(i)] = calculate_Q(par_IAV_ini, par_consider_ind, par_consider(i, :), y_ini);
end

delete(p);
delete(gcp('nocreate'))

% write Qs to files
file_QV = fopen('QV.txt', 'w');
for i = 1:par_num
    fprintf(file_QV, '%f\n', QVs(i));
end
fclose(file_QV);

file_QM = fopen('QM.txt', 'w');
for i = 1:par_num
    fprintf(file_QM, '%f\n', QMs(i));
end
fclose(file_QM);

file_QK = fopen('QK.txt', 'w');
for i = 1:par_num
    fprintf(file_QK, '%f\n', QKs(i));
end
fclose(file_QK);

toc;
disp(['run time:', num2str(toc)]);
