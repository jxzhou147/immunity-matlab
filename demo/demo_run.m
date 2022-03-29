clear;clc

figure;
for theta = -pi/2:pi:pi/2
    v0 = 1000;
    e0 = 1e6;

    par = [2e-2; e0; 10 ^ 6; 10; 1e-6; 1; 2.7e3; theta];

    tspan = 0:0.1:10;
    y0 = [v0; e0];
    [t, y] = ode45(@demo, tspan, y0, [], par);

%     plot(t, y(:, 1));
    plot(t, y(:, 1), 'LineWidth', 2);
    set(gca,'Fontsize',26); box on;
    xlabel('t (day)'); ylabel('V');
    hold on;
    pause(1);
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
    
    dydt = [0; 0];
    
    v = y(1);
    e = y(2);
    
    ce = par(1);
    e0 = par(2);
    kv = par(3);
    p = par(4);
    cv = par(5);
    r = par(6);
    ke = par(7);
    
    theta = par(8);
    eta = 1 * (sin(2*pi * t + theta)) +1;
    
%     if (t > 0.25)
%         eta = 1;
%     end
%     
%     
    if v < 1e-1
        v = 0;
        dydt(1) = 0;
    else
%         dydt(1) = (1 + (theta - 1) * 0.4) * p * v * (1 - v / kv) - cv * v * e;
        dydt(1) =  eta * p * v * (1 - v / kv) - cv * v * e;
    end
    
    if e < 1e-1
        e = 0;
        dydt(2) = 0;
    else
        dydt(2) = eta * r * e * (v / (v + ke)) - ce * e + ce * e0;
    end
    
end