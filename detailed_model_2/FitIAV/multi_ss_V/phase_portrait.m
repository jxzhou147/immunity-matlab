% calculate and plot zoros-root lines
par_set = importdata('multi_ss_par.txt');
par_set = par_set(50, :);

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
[CCL2, Mono] = meshgrid(-1e3:1e1:1e3, -1e3:1e1:1e3);
hadl = streamslice(CCL2, Mono, dMonodt(Mono, CCL2), dCCL2dt(Mono, CCL2));
xlabel('CCL2'); ylabel('Mono');
set(gca, 'XLim', [-1e3 1e3], 'YLim', [-1e3 1e3], 'Fontsize', 26, 'linewidth', 2);
hold on;

% plot zero-solution lines
CCL2 = (-1e3:10:1e3);
plot(CCL2, Mono1(CCL2), 'r', CCL2, Mono2(CCL2), 'g');
hold off;