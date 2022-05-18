clear;clc

% load Sis
S1s_V = load('S1s_V.txt');
S1s_M = load('S1s_M.txt');
S1s_K = load('S1s_K.txt');

STs_V = load('STs_V.txt');
STs_M = load('STs_M.txt');
STs_K = load('STs_K.txt');

% plot the time courses of Sobol indices
par_num = 7;
bin_num = 50;

rgb_1V = zeros(par_num, bin_num, 3);
rgb_1M = zeros(par_num, bin_num, 3);
rgb_1K = zeros(par_num, bin_num, 3);
rgb_TV = zeros(par_num, bin_num, 3);
rgb_TM = zeros(par_num, bin_num, 3);
rgb_TK = zeros(par_num, bin_num, 3);

for i = 1:par_num
    [rgb_1V(i, :, 1), rgb_1V(i, :, 2), rgb_1V(i, :, 3)] = ...
        arr2rgb(S1s_V(:, i), bin_num, 3, 116, 189);
    [rgb_1M(i, :, 1), rgb_1M(i, :, 2), rgb_1M(i, :, 3)] = ...
        arr2rgb(S1s_M(:, i), bin_num, 3, 116, 189);
    [rgb_1K(i, :, 1), rgb_1K(i, :, 2), rgb_1K(i, :, 3)] = ...
        arr2rgb(S1s_K(:, i), bin_num, 3, 116, 189);
    
    [rgb_TV(i, :, 1), rgb_TV(i, :, 2), rgb_TV(i, :, 3)] = ...
        arr2rgb(STs_V(:, i), bin_num, 216, 82, 24);
    [rgb_TM(i, :, 1), rgb_TM(i, :, 2), rgb_TM(i, :, 3)] = ...
        arr2rgb(STs_M(:, i), bin_num, 216, 82, 24);
    [rgb_TK(i, :, 1), rgb_TK(i, :, 2), rgb_TK(i, :, 3)] = ...
        arr2rgb(STs_K(:, i), bin_num, 216, 82, 24);
end

figure;
hold on
for i = 1:par_num
    bar_data_1V = ones(1, bin_num) * ...
        (max(S1s_V(:, i)) - min(S1s_V(:, i))) / bin_num;
    b = barh([4*i-1, 4*i], [bar_data_1V; zeros(1, bin_num)], 0.4, 'stacked');
    for j = 1:bin_num
        set(b(j), 'FaceColor', rgb_1V(i, j, :), 'EdgeColor', 'none');
    end
    
    bar_data_TV = ones(1, bin_num) * ...
        (max(STs_V(:, i)) - min(STs_V(:, i))) / bin_num;
    b = barh([4*i, 4*i+1], [bar_data_TV; zeros(1, bin_num)], 0.4, 'stacked');
    for j = 1:bin_num
        set(b(j), 'FaceColor', rgb_TV(i, j, :), 'EdgeColor', 'none');
    end
end
box on
set(gca, 'XLim', [-0.1 1.1], 'YLim', [0 31]);
xticks([0 0.5 1]);
yticks([3.5, 7.5, 11.5, 15.5, 19.5, 23.5, 27.5]);
yticklabels({'a_{MI}', 'a_{KI}', '\gamma', 'c_{M_IL6}', 'c_{M_IL10}', 'c_{M_CCL2}', 'c_{I_CXCL5}'});
xlabel('Sobol index');
title('V');
hold off;

figure;
hold on
for i = 1:par_num
    bar_data_1M = ones(1, bin_num) * ...
        (max(S1s_M(:, i)) - min(S1s_M(:, i))) / bin_num;
    b = barh([4*i-1, 4*i], [bar_data_1M; zeros(1, bin_num)], 0.4, 'stacked');
    for j = 1:bin_num
        set(b(j), 'FaceColor', rgb_1M(i, j, :), 'EdgeColor', 'none');
    end
    
    bar_data_TM = ones(1, bin_num) * ...
        (max(STs_M(:, i)) - min(STs_M(:, i))) / bin_num;
    b = barh([4*i, 4*i+1], [bar_data_TM; zeros(1, bin_num)], 0.4, 'stacked');
    for j = 1:bin_num
        set(b(j), 'FaceColor', rgb_TM(i, j, :), 'EdgeColor', 'none');
    end
end
box on
set(gca, 'XLim', [-0.1 1.1], 'YLim', [0 31]);
xticks([0 0.5 1]);
yticks([3.5, 7.5, 11.5, 15.5, 19.5, 23.5, 27.5]);
yticklabels({'a_{MI}', 'a_{KI}', '\gamma', 'c_{M_IL6}', 'c_{M_IL10}', 'c_{M_CCL2}', 'c_{I_CXCL5}'});
xlabel('Sobol index');
title('M');
hold off;

figure;
hold on
for i = 1:par_num
    bar_data_1K = ones(1, bin_num) * ...
        (max(S1s_K(:, i)) - min(S1s_K(:, i))) / bin_num;
    b = barh([4*i-1, 4*i], [bar_data_1K; zeros(1, bin_num)], 0.4, 'stacked');
    for j = 1:bin_num
        set(b(j), 'FaceColor', rgb_1K(i, j, :), 'EdgeColor', 'none');
    end
    
    bar_data_TK = ones(1, bin_num) * ...
        (max(STs_K(:, i)) - min(STs_K(:, i))) / bin_num;
    b = barh([4*i, 4*i+1], [bar_data_TK; zeros(1, bin_num)], 0.4, 'stacked');
    for j = 1:bin_num
        set(b(j), 'FaceColor', rgb_TK(i, j, :), 'EdgeColor', 'none');
    end
end
box on
set(gca, 'XLim', [-0.1 1.1], 'YLim', [0 31]);
xticks([0 0.5 1]);
yticks([3.5, 7.5, 11.5, 15.5, 19.5, 23.5, 27.5]);
yticklabels({'a_{MI}', 'a_{KI}', '\gamma', 'c_{M_IL6}', 'c_{M_IL10}', 'c_{M_CCL2}', 'c_{I_CXCL5}'});
xlabel('Sobol index');
title('K');
hold off;