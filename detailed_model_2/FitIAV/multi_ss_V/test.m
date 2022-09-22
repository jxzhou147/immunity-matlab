% load fitted parameter set ZT11
par_base = importdata('lhs_par.txt');
par_base = par_base(6, :);

% translate parFit to par in the odes
% log_par_ind = [1:39 43:54];
% for i = log_par_ind
%     par_base(i) = 10 .^ par_base(i);
% end

% solve odes
tmax = 50000;
tspan = 0:1:tmax;
y0 = zeros(14, 1);
y0(3) = 0;
y0(1) = par_base(55);
y0(4) = par_base(56);
y0(7) = par_base(57);
y0(12) = par_base(58);
y0(13) = 10;

[t, y] = ode15s(@ODE_IAV, tspan, y0, [], par_base);
y = real(y);