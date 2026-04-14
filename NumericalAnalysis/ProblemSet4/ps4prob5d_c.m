clear; clc;

yexact = @(t) 0.5*(exp(t) - sin(t) - cos(t));

kmax = 16;                  % gives Nh = 10,20,40,...,10*2^15
Nh_vals = 10 * 2.^(0:kmax-1);

h_vals = 1 ./ Nh_vals;
errors = zeros(size(Nh_vals));

for k = 1:length(Nh_vals)
    Nh = Nh_vals(k);
    h = 1 / Nh;

    t = 0;
    y = 0;

    for n = 1:Nh
        tnew = t + h;
        y = (y + h*sin(tnew)) / (1 - h);
        t = tnew;
    end

    errors(k) = abs(y - yexact(1));
end

loglog(h_vals, errors, 'o-','LineWidth',2);
grid on;
xlabel('h');
ylabel('Error at t=1');
title('Backward Euler Error vs Step Size');

hold on;
loglog(h_vals, h_vals, '--');
legend('Backward Euler error','O(h) reference','Location','best');

