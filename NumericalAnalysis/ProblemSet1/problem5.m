x = linspace(-5e-6,5e-6,1001);
y = log(1+x) - cos(x) - x +1;
z = (log(1+x) - x) + (1 - cos(x));
plot(x,y,'.k')
hold on
plot(x,z,'.b')
legend('log(1+x) - cos(x) - x +1', '(log(1+x) - x) + (1 - cos(x))')
hold off