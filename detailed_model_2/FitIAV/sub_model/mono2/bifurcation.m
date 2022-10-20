% plot bifurcation diagram

% fix 1 parameter set and let I be the consider par
% par_base = importdata('multi_ss_par.txt');
% par_base = par_base(8, :);

% fix 1 parameter set and let c_I_CCL2 be the consider par
par_base = importdata('par_M_init.txt');
par_base = par_base.data';
% translate to par in the odes
% log_par_ind = (1:7);
% for i = log_par_ind
%     par_base(i) = 10 .^ par_base(i);
% end

par_set_num = 100;
par_consider_idx = 6;
par_consider = zeros(par_set_num, 1);
multi_ss_bool = false(par_set_num, 1);
multi_ss = zeros(2, 2, par_set_num);

p = parpool(20);
tic;
for i = 1:par_set_num
    par_consider(i) = i * 5e-2;
    [multi_ss_bool(i), multi_ss(:, :, i)] = if_multi_ss(par_base, par_consider_idx, par_consider(i));
end
delete(p);
toc;

disp(['run time:', num2str(toc)]);
disp(['number of multi_ss: ', num2str(sum(multi_ss_bool))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot figures
figure;
xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);
hold on;
subplot(1,2,1); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, reshape(multi_ss(1, 1, :), 100, 1), 'r', 'LineWidth', 2);
plot(par_consider, reshape(multi_ss(2, 1, :), 100, 1), 'b', 'LineWidth', 2);
xlabel('c_{I-CCL2}'); ylabel('Mono');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold on;

subplot(1,2,2); hold on; set(gca,'Fontsize',26); box on;
plot(par_consider, reshape(multi_ss(1, 2, :), 100, 1), 'r', 'LineWidth', 2);
plot(par_consider, reshape(multi_ss(2, 2, :), 100, 1), 'b', 'LineWidth', 2);
xlabel('c_{I-CCL2}'); ylabel('CCL2');
% set(gca, 'XTick', [100:150:400], 'XLim', [100 400], 'YLim', [-0.6 8], 'Fontsize', 26, 'linewidth', 2);
hold off;