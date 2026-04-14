clear; clc;

f = @(t,y) sin(t) + y;
yexact = @(t) 0.5*(exp(t) - sin(t) - cos(t));

kmax = 20;                 % up to Nh = 2^20
Nh_vals = 2.^(1:kmax);     % Nh = 2,4,8,...,2^20

h_vals = 1 ./ Nh_vals;     % step sizes
errors = zeros(size(Nh_vals));

for k = 1:length(Nh_vals)
    Nh = Nh_vals(k);
    h = 1 / Nh;

    t = 0;
    y = 0;

    % Forward Euler loop
    for n = 1:Nh
        y = y + h*(sin(t) + y);
        t = t + h;
    end

    % compute error at t=1
    errors(k) = abs(y - yexact(1));
end

% log-log plot
loglog(h_vals, errors, 'o-','LineWidth',2);
grid on;
xlabel('h');
ylabel('Error at t=1');
title('Forward Euler Error vs Step Size (log-log)');

% optional: reference slope line ~ h
hold on;
loglog(h_vals, h_vals, '--'); % slope 1 reference
legend('Error','O(h) reference','Location','best');