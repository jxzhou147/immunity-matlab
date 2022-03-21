clear;clc

for eta = 0.5:0.1:0.7
    v0 = 10;
    e0 = 1e6;
    m0 = 100;
    %eta = 1;

    par = [2e-2; 2e-2 * e0; 10 ^ 5.5; 4.4; 1.24e-6; 0.33; 2.7e3; eta];

    tspan = 0:0.1:10;
    y0 = [v0; e0; m0];
    [t, y] = ode45(@demo, tspan, y0, [], par);

%     plot(y(:, 2), y(:, 1));
    plot(t, log10(y(:, 1)));
    hold on;
    pause(0.1);
    drawnow;
end

% figure;
% xSize = 20; X=xSize; ySize = 7;xLeft = (xSize-xSize)/2; Y=ySize; yTop = (ySize-ySize)/2;
% set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);set(gcf,'Position',[X Y xSize*50 ySize*55]);

% subplot(1,2,1);
% plot(t, log10(y(:, 1)), 'r', 'LineWidth', 2);
% subplot(1,2,2);
% plot(t, y(:, 2), 'b', 'LineWidth', 2);




function dydt = demo(t, y, par)
    
    dydt = [0; 0; 0];
    
    v = y(1);
    e = y(2);
    m = y(3);
    
    ce = par(1);
    se = par(2);
    kv = par(3);
    p = par(4);
    cv = par(5);
    r = par(6);
    ke = par(7);
    eta = par(8);
    
%     if (t > 3)
%         eta = 1;
%     end
    
    q = 1e2;
    km = 100;
    cm = 1;
    sm = 0;
    gamma = 1e5;
    km1 = 1e5;
    
    if v < 1e-1
        v = 0;
        dydt(1) = 0;
    else
        dydt(1) = p * v * (1 - v / kv) - cv * v * e;
    end
    
    if e < 1e-1
        e = 0;
        dydt(2) = 0;
    else
        dydt(2) = r * e * (v / (v + ke)) - ce * e + se;
    end
    
    if m < 1e-1
        m = 0;
        dydt(3) = 0;
    else
        dydt(3) = eta * q * (v / (v + km)) + gamma * m / (m + km1) - cm * m + cm * 100;
    end
    
end