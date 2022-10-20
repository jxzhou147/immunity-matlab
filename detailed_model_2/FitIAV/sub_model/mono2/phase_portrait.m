clear;clc

% calculate and plot zoros-root lines
% par_set = importdata('multi_ss_par.txt');
% par_set = par_set(1, :);
par_set = importdata('par_M_init.txt');
par_set = par_set.data';
% translate to par in the odes
log_par_ind = (1:7);
for i = log_par_ind
    par_set(i) = 10 .^ par_set(i);
end

par_set = num2cell(par_set, 1);
[c_CCL2_Mono, K_CCL2_Mono, n_CCL2_Mono, d_Mono, c_M_CCL2, c_I_CCL2, d_CCL2, I] = deal(par_set{:});

% ODEs
dMonodt = @(Mono, CCL2)c_CCL2_Mono .* CCL2 .^ n_CCL2_Mono ./ ...
    (K_CCL2_Mono .^ n_CCL2_Mono + CCL2 .^ n_CCL2_Mono) - d_Mono .* Mono;
    
dCCL2dt = @(Mono, CCL2)c_M_CCL2 .* Mono + c_I_CCL2 .* I - d_CCL2 .* CCL2;

% zero-solution lines
Mono1 = @(CCL2)(c_CCL2_Mono .* CCL2 .^ n_CCL2_Mono ./ ...
    (K_CCL2_Mono .^ n_CCL2_Mono + CCL2 .^ n_CCL2_Mono) ./ d_Mono);

Mono2 = @(CCL2)((d_CCL2 .* CCL2 - c_I_CCL2 .* I) ./ c_M_CCL2);

% plot phase portrait
% [CCL2, Mono] = meshgrid(0:1e1:1500, 0:1e1:400);
% hadl = streamslice(CCL2, Mono, dMonodt(Mono, CCL2), dCCL2dt(Mono, CCL2));
% xlabel('CCL2'); ylabel('Mono');
% set(gca, 'XLim', [0 1500], 'YLim', [0 400], 'Fontsize', 26, 'linewidth', 2);
% hold on;

% plot zero-solution lines
CCL2 = (0:10:1000);
plot(CCL2, Mono1(CCL2), 'r', 'LineWidth', 2); hold on;
plot(CCL2, Mono2(CCL2), 'b', 'LineWidth', 2); hold on;
set(gca, 'LineWidth', 2, 'FontSize', 26); hold on;
xlabel('CCL2'); ylabel('Mono');
hold off;