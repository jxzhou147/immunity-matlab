clear;clc

% load clock model parameters
par_clock = load('par_clock.csv');
dp_rev = par_clock(8);

tspan = 0:1:500;
y0 = zeros(12, 1);

% Loop for checking dp_rev's effects on period and amplitude
num = 41;
amp = zeros(num, 13);
amp_fold = zeros(num, 13); % amp / (amp of ratio = 1)
period = zeros(num, 13);
for i = 1:num
    ratio = 1 + (i - 1) * 0.1;
    par_clock(8) = dp_rev * ratio;
    
    % solve odes
    [t, y] = ode15s(@ODE_Clock, tspan, y0, [], par_clock);
    y = real(y);
    
    for j = 1:12
        [pks, locs] = findpeaks(y(200:500, j), tspan(200:500));
        [troughs, locs_tr] = findpeaks(-1 * y(200:500, j), tspan(200:500));
        troughs =  -1 * troughs;
        amp(i, j) = mean(pks) - mean(troughs);
        amp_fold(i, j) = amp(i, j) ./ amp(1, j);
        period(i, j) = mean(diff(locs));
    end
    amp(i, 13) = ratio;
    amp_fold(i, 13) = ratio;
    period(i, 13) = ratio;
end

% figures
figure;
plot(period(:, 13), period(:, 1:12), 'LineWidth', 2); hold on;
xlabel('fold change of dp_{rev}'); ylabel('period (h)');
set(gca, 'XLim', [1 5], 'Fontsize', 26);
hold off;

figure;
plot(amp(:, 13), amp(:, 1:12), 'LineWidth', 2); hold on;
xlabel('fold change of dp_{rev}'); ylabel('amplitude');
set(gca, 'XLim', [1 5], 'Fontsize', 26);
hold off;

figure;
plot(amp_fold(:, 13), amp_fold(:, 1:12), 'LineWidth', 2); hold on;
xlabel('fold change of dp_{rev}'); ylabel('fold change of amplitude');
set(gca, 'XLim', [1 5], 'Fontsize', 26);
hold off;