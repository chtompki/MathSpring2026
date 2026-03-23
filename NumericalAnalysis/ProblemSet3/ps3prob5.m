% finite difference approximations for derivatives
format long
clc; clear;

%% part (b): centered first derivative for exp(x) at x=0
f = @(x) exp(x);
fp1 = @(x) exp(x);

h = 10.^(-1:-1:-5);
x0 = 0.0;
fp1_exct = fp1(x0);

fp1_cnfd2 = zeros(size(h));
err_fp1_cnfd2 = zeros(size(h));

for i = 1:length(h)
    hi = h(i);
    c = fdcoeffV(1,[-hi hi]);
    fp1_cnfd2(i) = c(1)*f(x0-hi) + c(2)*f(x0+hi);
    err_fp1_cnfd2(i) = abs(fp1_cnfd2(i)-fp1_exct);
end

figure;
loglog(h,err_fp1_cnfd2,'o');
xlabel('h')
ylabel('error')
title('Centered first derivative, f(x)=e^x at x=0')

c = polyfit(log(h),log(err_fp1_cnfd2),1);
conv_order_fp1_cnfd2 = c(1)

%% part (c): fourth-order forward first derivative for sin(x) at x=0
f = @(x) sin(x);
fp1 = @(x) cos(x);

h = 10.^(-1:-1:-5);
x0 = 0.0;
fp1_exct = fp1(x0);

fp1_fwfd4 = zeros(size(h));
err_fp1_fwfd4 = zeros(size(h));

for i = 1:length(h)
    hi = h(i);
    c = fdcoeffV(1,[0 hi 2*hi 3*hi 4*hi]);
    fp1_fwfd4(i) = c(1)*f(x0) ...
                 + c(2)*f(x0+hi) ...
                 + c(3)*f(x0+2*hi) ...
                 + c(4)*f(x0+3*hi) ...
                 + c(5)*f(x0+4*hi);
    err_fp1_fwfd4(i) = abs(fp1_fwfd4(i)-fp1_exct);
end

figure;
loglog(h,err_fp1_fwfd4,'s');
xlabel('h')
ylabel('error')
title('Fourth-order forward first derivative, f(x)=sin(x) at x=0')

c = polyfit(log(h),log(err_fp1_fwfd4),1);
conv_order_fp1_fwfd4 = c(1)