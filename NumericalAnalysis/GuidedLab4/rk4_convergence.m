clear;
clc;

f = @(t,y) 1 - y + exp(2*t).*y.^2;
exact = @(t) exp(-t).*tan(exp(t)-1);

a = 0;
b = 0.9;
alpha = 0;

% Moderate step sizes using the original rk4.m.
hvals = [0.1 0.05 0.025 0.0125 0.00625 0.003125 0.0015625];
errors = zeros(size(hvals));

for j = 1:length(hvals)
    h = hvals(j);
    N = round((b-a)/h);

    [~,w] = rk4(f,a,b,alpha,N);
    errors(j) = abs(exact(b)-w(end));
end

disp('       h                 error')
disp([hvals(:),errors(:)])

% Estimate the slope using the last four values, where the method is
% close to its asymptotic fourth-order regime.
p = polyfit(log(hvals(end-3:end)),log(errors(end-3:end)),1);
fprintf('Observed convergence rate = %.6f\n',p(1));

figure;
loglog(hvals,errors,'o-','LineWidth',1.5);
grid on;
xlabel('Stepsize h');
ylabel('Error at t = 0.9');
title('RK4 Convergence');
saveas(gcf,'rk4_convergence.png');
