f = @(t,y) [y(2); -y(1)];

alpha = [0;1];

[t,w] = rk4(f,0,2*pi,alpha,100);

plot(t,w(:,1),'LineWidth',1.5)
hold on
plot(t,w(:,2),'LineWidth',1.5)
legend('y_1','y_2')
xlabel('t')
ylabel('Approximate solution')
grid on



f = @(t,y) y - t.^2 + 1;

[t,w] = rk4(f,0,2,0.5,10);

disp([t,w])