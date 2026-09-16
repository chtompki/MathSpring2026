%% Problem 4: absolute stability of classical RK4.
% Run in MATLAB; no toolboxes or cluster required.
% All figures and data are saved beside this script.
clearvars; clc;
outdir = fileparts(mfilename('fullpath'));
Q = @(z) 1+z+z.^2/2+z.^3/6+z.^4/24;
T = 10;
beta = fzero(@(r) Q(-r)-1,[2,3]);
lambda = -55.31;
hcrit = beta/abs(lambda);
fprintf('Negative-real boundary = %.12f; hcrit = %.12f\n',-beta,hcrit);

%% Parts (b) and (d): contour |Q(z)|=1 and shaded |Q(z)|<=1.
x = linspace(-3.5,1.5,1001);
y = linspace(-3.5,3.5,1001);
[X,Y] = meshgrid(x,y); Z = X+1i*Y;
A = abs(Q(Z));
fig = figure('Color','w','Visible','off');
[C,boundary] = contour(X,Y,A,[1,1],'k','LineWidth',1.5); hold on;
k=1;
while k<size(C,2)
    n=C(2,k); vertices=C(:,k+(1:n));
    patch(vertices(1,:),vertices(2,:),[0.83,0.90,0.98], ...
        'EdgeColor','none','HandleVisibility','off');
    k=k+n+1;
end
uistack(boundary,'top');
xline(0,'k:'); yline(0,'k:');
zComplex = 0.1*(1+20i);
point = plot(real(zComplex),imag(zComplex),'ro', ...
    'MarkerFaceColor','r','MarkerSize',7);
axis equal; xlim([-3.5,1.5]); ylim([-3.5,3.5]); grid on;
xlabel('Re(z), z=h lambda'); ylabel('Im(z)');
title('Classical RK4 absolute stability region');
legend([boundary,point],{'|Q(z)|=1','z=0.1+2i'},'Location','southwest');
exportgraphics(fig,fullfile(outdir,'stability_region.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outdir,'stability_region.png'),'Resolution',180);
% Export each contour component for a self-contained LaTeX plot.
k=1; component=0;
while k<size(C,2)
    n=C(2,k); component=component+1;
    writematrix(C(:,k+(1:n))',fullfile(outdir,sprintf('boundary_%d.csv',component)));
    k=k+n+1;
end

%% Part (c): h=0.05 is stable and h=0.06 is unstable.
% h=0.06 does not divide 10. Use 166 full steps to 9.96, then 0.04.
hvalues=[0.05,0.06]; realRuns=struct([]);
fig = figure('Color','w','Visible','off');
for j=1:2
    h=hvalues(j);
    [t,w,steps]=rk4test(lambda,T,h);
    exact=exp(lambda*t);
    % Product of the RK4 amplification factors, including the last step.
    predicted=[1;cumprod(Q(lambda*steps))];
    assert(max(abs(w-predicted)./max(1,abs(predicted)))<1e-11);
    realRuns(j).h=h; realRuns(j).t=t; realRuns(j).w=w;
    realRuns(j).steps=steps; realRuns(j).exact=exact;
    fprintf('h=%.2f: Q=%.12f, w(10)=%.12e, last step=%.5f\n', ...
        h,Q(lambda*h),w(end),steps(end));
    semilogy(t,abs(w),'LineWidth',1.5); hold on;
    writetable(table(t,w,exact),fullfile(outdir,sprintf('real_h%03d.csv',round(1000*h))));
end
tExact=linspace(0,T,2001)';
semilogy(tExact,exp(lambda*tExact),'k--','LineWidth',1.2);
grid on; xlabel('t'); ylabel('Magnitude');
title('lambda=-55.31: stable and unstable RK4 steps');
legend('RK4 h=0.05','RK4 h=0.06 (last step 0.04)','Exact','Location','southwest');
exportgraphics(fig,fullfile(outdir,'negative_real_test.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outdir,'negative_real_test.png'),'Resolution',180);

%% Part (d): stable numerical decay can contradict physical growth.
lambdaComplex=1+20i; hComplex=0.1;
[t,w,steps]=rk4test(lambdaComplex,T,hComplex);
exact=exp(lambdaComplex*t);
predicted=[1;cumprod(Q(lambdaComplex*steps))];
assert(max(abs(w-predicted))<1e-12);
q=Q(zComplex);
fprintf('Q(0.1+2i)=%.12f%+.12fi; magnitude=%.12f\n',real(q),imag(q),abs(q));
fprintf('At t=10: |w|=%.12e; |exact|=%.12e\n',abs(w(end)),abs(exact(end)));
complexRun=struct('t',t,'w',w,'exact',exact,'Q',q);
writetable(table(t,real(w),imag(w),abs(w),abs(exact), ...
    'VariableNames',{'t','real_w','imag_w','abs_w','abs_exact'}), ...
    fullfile(outdir,'complex_test.csv'));
fig=figure('Color','w','Visible','off');
semilogy(t,abs(w),'b-','LineWidth',1.5); hold on;
semilogy(t,abs(exact),'k--','LineWidth',1.5);
grid on; xlabel('t'); ylabel('Magnitude');
title('lambda=1+20i, h=0.1: false numerical damping');
legend('|w_n|','|exp(lambda t)|=exp(t)','Location','best');
exportgraphics(fig,fullfile(outdir,'complex_growth_test.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outdir,'complex_growth_test.png'),'Resolution',180);
save(fullfile(outdir,'stability_results.mat'), ...
    'beta','hcrit','realRuns','complexRun');

function [t,w,steps]=rk4test(lambda,T,h)
    N=floor(T/h);
    t=(0:N)'*h;
    if T-t(end)>100*eps(max(1,T))
        t=[t;T];
    else
        t(end)=T;
    end
    steps=diff(t);
    w=complex(zeros(size(t))); w(1)=1;
    % Actual four-stage RK4, not just multiplication by Q.
    for n=1:numel(steps)
        hn=steps(n);
        k1=lambda*w(n);
        k2=lambda*(w(n)+hn*k1/2);
        k3=lambda*(w(n)+hn*k2/2);
        k4=lambda*(w(n)+hn*k3);
        w(n+1)=w(n)+hn*(k1+2*k2+2*k3+k4)/6;
    end
    if isreal(lambda), w=real(w); end
end

% Negative-real boundary = -2.785293563405; hcrit = 0.050357865909
% h=0.05: Q=0.970565404941, w(10)=2.540747063848e-03, last step=0.05000
% h=0.06: Q=2.150291693782, w(10)=6.704798852280e+54, last step=0.04000
% Q(0.1+2i)=-0.438162500000+0.743666666667i; magnitude=0.863149168752
% At t=10: |w|=4.060557481727e-07; |exact|=2.202646579481e+04