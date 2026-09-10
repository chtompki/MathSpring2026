clear;
clc;

f = @(t,y) 1 - y + exp(2*t).*y.^2;
exact = @(t) exp(-t).*tan(exp(t)-1);

a = 0;
b = 0.9;
alpha = 0;
h = 0.1;
N = round((b-a)/h);

[t,w] = rk4(f,a,b,alpha,N);

x = exact(t);
err = x - w;

disp('        t              x              w             err')
disp([t,x,w,err])

figure;
plot(t,x,'-','LineWidth',1.5);
hold on;
plot(t,w,'o--','LineWidth',1.5);
grid on;
xlabel('t');
ylabel('y');
legend('Exact solution','RK4 approximation','Location','best');
title('Exact Solution and RK4 Approximation, h = 0.1');
