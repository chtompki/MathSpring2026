%% Problem 3: reproduce Figure 7.21 using classical RK4.
% Run this file in MATLAB. No cluster or toolboxes are needed.
% Outputs are saved beside this script; results remain in the workspace.
clearvars; clc;
outdir = fileparts(mfilename('fullpath'));
m = 16; a = 0.25; T = 12; y0 = [0.5;0.02];
hvalues = [0.05,0.06];
F = @(t,y) [1-y(1)-m*y(1)*y(2)/(a+y(1)); ...
                 m*y(1)*y(2)/(a+y(1))-y(2)];
E0 = [1;0];
Estar = [a/(m-1);1-a/(m-1)];
J = @(s,x) [-1-m*a*x/(a+s)^2, -m*s/(a+s); ...
                 m*a*x/(a+s)^2, m*s/(a+s)-1];
lambda0 = eig(J(E0(1),E0(2)));
lambdaStar = eig(J(Estar(1),Estar(2)));
R = @(z) 1+z+z.^2/2+z.^3/6+z.^4/24;
beta = fzero(@(b) R(-b)-1,[2,3]);
hcrit = beta/max(abs(lambdaStar));
fprintf('Washout: (%g,%g), eigenvalues %g and %g\n',E0,lambda0);
fprintf('Coexistence: (%.12f,%.12f), eigenvalues %g and %g\n', ...
    Estar,lambdaStar);
fprintf('RK4 local stability limit: h < %.12f\n',hcrit);

fig = figure('Color','w','Position',[100,100,960,800]);
tiledlayout(2,1,'TileSpacing','compact','Padding','compact');
results = struct([]);
summaryData = zeros(2,8);
for j = 1:numel(hvalues)
    h = hvalues(j);
    [t,W] = rk4system(F,T,y0,h);
    s = W(:,1); x = W(:,2);
    % Exact total mass follows q'=1-q, q(0)=0.52.
    qexact = 1+(sum(y0)-1)*exp(-t);
    massError = max(abs(s+x-qexact));
    % RK4 also advances this linear mass balance with factor R(-h).
    qRK4 = 1+(sum(y0)-1)*R(-h).^(0:numel(t)-1)';
    massIdentityError = max(abs(s+x-qRK4));
    assert(all(isfinite(W(:))),'Nonfinite RK4 values.');
    assert(massIdentityError<1e-11,'Discrete mass-balance check failed.');
    tail = t>=10;
    amp = R(h*lambdaStar);
    summaryData(j,:) = [h,s(end),x(end),min(s(tail)),max(s(tail)), ...
        massError,max(abs(amp)),massIdentityError];
    results(j).h = h; results(j).t = t; results(j).w = W;
    results(j).qexact = qexact; results(j).amplification = amp;
    results(j).massError = massError;
    data = table(t,s,x,qexact,s+x-qexact, ...
        'VariableNames',{'t','s','x','q_exact','mass_error'});
    writetable(data,fullfile(outdir,sprintf('chemostat_h%03d.csv',round(1000*h))));
    fprintf('\nh=%.2f: s(12)=%.12f, x(12)=%.12f\n',h,s(end),x(end));
    fprintf('  RK4 eigenmode amplification factors: %.12f, %.12f\n',amp);
    fprintf('  Maximum total-mass error: %.3e\n',massError);
    fprintf('  Late substrate range (10<=t<=12): [%.12f, %.12f]\n', ...
        min(s(tail)),max(s(tail)));
    nexttile;
    plot(t,s,'k-','LineWidth',1.2); hold on;
    plot(t,x,'k:','LineWidth',1.6);
    xlim([0,T]); ylim([0,1.05]); xticks(0:2:T); yticks(0:0.2:1);
    xlabel('Nondimensional time'); ylabel('Nondimensional mass');
    title(sprintf('RK4 with h = %.2f',h));
    legend('Substrate s(t)','Cell mass x(t)','Location','east');
    box on;
end
summary = array2table(summaryData,'VariableNames', ...
    {'h','s_final','x_final','s_tail_min','s_tail_max', ...
     'max_mass_error','max_amplification','mass_identity_error'});
disp(summary);
writetable(summary,fullfile(outdir,'chemostat_summary.csv'));
exportgraphics(fig,fullfile(outdir,'chemostat_comparison.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outdir,'chemostat_comparison.png'),'Resolution',200);
save(fullfile(outdir,'chemostat_results.mat'), ...
    'results','summary','E0','Estar','lambda0','lambdaStar','hcrit');

%% A finer-step check of the h=0.05 trajectory at its grid points.
[tFine,WFine] = rk4system(F,T,y0,0.005);
coarseFineDifference = max(abs(results(1).w-WFine(1:10:end,:)),[],1);
fprintf('\nMax difference from h=0.005 on h=0.05 nodes:\n');
fprintf('  substrate %.6e; cell mass %.6e\n',coarseFineDifference);

function [t,W] = rk4system(F,T,y0,h)
    N = round(T/h);
    assert(abs(N*h-T)<100*eps(max(1,T)),'h must divide T.');
    t = (0:N)'*h;
    W = zeros(N+1,numel(y0));
    W(1,:) = y0(:)';
    for n = 1:N
        y = W(n,:)';
        k1 = F(t(n),y);
        k2 = F(t(n)+h/2,y+h*k1/2);
        k3 = F(t(n)+h/2,y+h*k2/2);
        k4 = F(t(n)+h,y+h*k3);
        W(n+1,:) = (y+h*(k1+2*k2+2*k3+k4)/6)';
    end
end


% Washout: (1,0), eigenvalues -1 and 11.8
% Coexistence: (0.016666666667,0.983333333333), eigenvalues -55.3125 and -1
% RK4 local stability limit: h < 0.050355589847
%
% h=0.05: s(12)=0.016666721334, x(12)=0.983330329442
%   RK4 eigenmode amplification factors: 0.970748744905, 0.951229427083
%   Maximum total-mass error: 9.589e-09
%   Late substrate range (10<=t<=12): [0.016666721334, 0.016667069114]
%
% h=0.06: s(12)=0.202874124682, x(12)=0.797122926092
%   RK4 eigenmode amplification factors: 2.150727239037, 0.941764540000
%   Maximum total-mass error: 2.005e-08
%   Late substrate range (10<=t<=12): [0.086861233814, 0.202874124682]
%      h      s_final     x_final    s_tail_min    s_tail_max    max_mass_error    max_amplification    mass_identity_error
%     ____    ________    _______    __________    __________    ______________    _________________    ___________________
%
%     0.05    0.016667    0.98333     0.016667      0.016667       9.5885e-09           0.97075             4.4409e-16
%     0.06     0.20287    0.79712     0.086861       0.20287       2.0045e-08            2.1507             2.2204e-16
%
%
% Max difference from h=0.005 on h=0.05 nodes:
%   substrate 2.272732e-03; cell mass 2.272740e-03