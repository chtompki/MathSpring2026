%% RK4 solution: y' = 1-y+exp(2*t)*y^2, y(0)=0, 0<=t<=0.9.
% Run: rk4_ivp_solution
% Default: full sweep h=10.^(-(1:9)), about one billion total steps.
% Quick check: setenv('RK4_QUICK','1'); rk4_ivp_solution
% Restore full mode: setenv('RK4_QUICK','0');
% Output files are written alongside this script. No toolboxes required.
clearvars; clc;
outdir = fileparts(mfilename('fullpath'));
T = 0.9;
exact = @(t) exp(-t).*tan(expm1(t));

%% Part 2: store every value, including the initial condition.
h = 0.1;
N = round(T/h);
t = (0:N)'*h;
w = zeros(N+1,1);
for n = 1:N
    w(n+1) = rk4step(t(n),w(n),h);
end
x = exact(t);
err = x-w;                         % signed errors, as requested
nodal = table(t,x,w,err);
disp(nodal);
writetable(nodal,fullfile(outdir,'rk4_nodes.csv'));
save(fullfile(outdir,'rk4_nodes.mat'),'t','x','w','err','h');

%% Parts 3 and challenge: constant-memory endpoint solver.
exponents = 1:9;
if strcmp(getenv('RK4_QUICK'),'1')
    exponents = 1:5;
end
h_requested = 10.^(-exponents(:));
Nsteps = round(T./h_requested);
h_values = T./Nsteps;               % endpoint reached in N steps
m = numel(h_values);
w_final = NaN(m,1); signed_error = NaN(m,1);
abs_error = NaN(m,1); seconds = NaN(m,1); rate = NaN(m,1);
x_final = exact(T);
for j = 1:m
    timer = tic;
    w_final(j) = rk4final(T,Nsteps(j));
    seconds(j) = toc(timer);
    signed_error(j) = x_final-w_final(j);
    abs_error(j) = abs(signed_error(j));
    if j>1 && abs_error(j)>0 && abs_error(j-1)>0
        rate(j) = log(abs_error(j-1)/abs_error(j)) / ...
                  log(h_values(j-1)/h_values(j));
    end
    fprintf('h=%.1e N=%.0f E=%.12e rate=%.6f seconds=%.2f\n', ...
        h_values(j),Nsteps(j),abs_error(j),rate(j),seconds(j));
    % Save completed runs immediately, even if a later run is interrupted.
    convergence = table(h_values(1:j),Nsteps(1:j),w_final(1:j), ...
        signed_error(1:j),abs_error(1:j),rate(1:j),seconds(1:j), ...
        'VariableNames',{'h','N','w_final','signed_error','abs_error','rate','seconds'});
    writetable(convergence,fullfile(outdir,'rk4_convergence.csv'));
    save(fullfile(outdir,'rk4_convergence.mat'),'convergence','x_final');
end

%% Plot measured absolute errors; zero errors cannot appear on log axes.
fig = figure('Visible','off','Color','w');
positive = abs_error>0;
loglog(h_values(positive),abs_error(positive),'o-','LineWidth',1.5);
hold on;
reference = abs_error(2)*(h_values/h_values(2)).^4;
loglog(h_values,reference,'k--','LineWidth',1.2);
grid on; xlabel('Step size h'); ylabel('|y(0.9)-w_N|');
title('Classical RK4: truncation error and round-off');
legend('Measured endpoint error','C h^4','Location','northwest');
exportgraphics(fig,fullfile(outdir,'rk4_convergence.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outdir,'rk4_convergence.png'),'Resolution',200);
close(fig);
if any(~positive)
    fprintf('Zero computed errors were omitted from the log-log plot.\n');
end

function y = rk4final(T,N)
    h = T/N;
    y = 0;
    % A scalar while-loop guarantees no O(N) index/time/solution arrays.
    n = 0;
    while n<N
        t = n*h;                   % avoid accumulating t=t+h round-off
        y = rk4step(t,y,h);
        n = n+1;
    end
end

function next = rk4step(t,y,h)
    % Midpoint exponential is shared by stages 2 and 3.
    a = exp(2*t); b = exp(2*(t+h/2)); c = exp(2*(t+h));
    k1 = 1-y+a*y*y;
    z = y+h*k1/2; k2 = 1-z+b*z*z;
    z = y+h*k2/2; k3 = 1-z+b*z*z;
    z = y+h*k3;   k4 = 1-z+c*z*z;
    next = y+h*(k1+2*k2+2*k3+k4)/6;
end
