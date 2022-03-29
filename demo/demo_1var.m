clear;clc

omega = 2 * pi / 24;
I0 = 10;
k = 1;
d = 0.03;
t01 = 6;
t02 = 18;

par1 = [omega; t01; I0; d; k];
par2 = [omega; t02; I0; d; k];

tspan = 0:1:100;
y0 = 0;

[t, y1] = ode45(@demo1, tspan, y0, [], par1);
[t, y2] = ode45(@demo1, tspan, y0, [], par2);

plot(t, y1, t - 12, y2);


function dydt = demo1(t, y, par)
    
    omega = par(1);
    t0 = par(2);
    I0 = par(3);
    d = par(4);
    k = par(5);
    
    I = t >= t0;
    I = I * I0;
    
    if t < (t0 + 6)
        dydt = (sin(omega * t) ./ 2 + 1) ...
            * k * I - d * y;
    else
        dydt = k * I - d * y;
    end
end

