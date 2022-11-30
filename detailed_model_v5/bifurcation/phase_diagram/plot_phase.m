clear;clc

curve050 = importdata("curve050.dat");
curve060 = importdata("curve060.dat");
curve800 = importdata("curve800.dat");
curve900 = importdata("curve900.dat");

figure;
plot(curve050(:, 1), curve050(:, 2), 'k', 'LineWidth', 1); hold on;
plot(curve060(:, 1), curve060(:, 2), 'k', 'LineWidth', 1);
plot(curve800(:, 1), curve800(:, 2), 'k', 'LineWidth', 1);
plot(curve900(:, 1), curve900(:, 2), 'k', 'LineWidth', 1);
hold off;