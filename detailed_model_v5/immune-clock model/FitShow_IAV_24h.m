clear;clc;

% load immunity and clock-immunity parameters
par_IAV = importdata('par_base_IAV_clock_values.txt');
par_IAV = par_IAV.data;

% load clock parameters
par_clock = load('par_clock.csv');

% load inflammation parameters
par_infla = importdata('par_infla.txt');
par_infla = par_infla.data;

y0 = zeros(27, 1);
y0_IAV = importdata('init_base.txt');
y0_IAV = y0_IAV.data;
for i = 13:26
    y0(i) = y0_IAV(i - 12);
end

tmax = 2000;
tspan = 0:1:tmax;

ys = zeros(tmax+1, 27, 23);
t_IAV = zeros(23);

parfor i = 1:23
    t_IAV(i) = 493 + (i - 1);
    [~, ys(:, :, i)] = ode15s(@ODE_Clock_IAV, tspan, y0, [], par_clock, par_IAV, t_IAV(i), par_infla);
end

%%

clock_name = ["Per", "Cry", "Rev", "Ror", "Bmal1", "PER", "CRY", "REV", "ROR", "BMAL1", "PER-CRY", "CLOCK-BMAL1"];
figure;
for j = 1:12
    subplot(3, 4, j);
    hold on;
    for i = 1:23
        plot(tspan, ys(:, j, i));
    end
    xlim([490 550]);
    xlabel('Time (h)');
    ylabel(clock_name(j));
    hold off;
end

var_name = ["H", "If", "V", "M_0", "M", "Mono", "N", "IL1b", "IL10", "CCL2", "CXCL5", "K", "T", "T_E"];
figure;
for j = 13:26
    subplot(4, 4, j - 12);
    hold on;
    for i = 1:23
        plot(tspan - t_IAV(i), ys(:, j, i));
    end
    xlim([0 200]);
    ylabel(var_name(j - 12)); hold off;
end
