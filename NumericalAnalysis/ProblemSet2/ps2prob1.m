determinants=zeros(1,100);
cpuTime=zeros(1,100);
for n=100:100:10000
    A = magic(n);
    tstart = cputime;
    det(A);
    cpuTime(n/100) = cputime - tstart;
end
plotN=1:100;
figure;hold on; semilogy(plotN, 0.000025.*plotN.^3, '--x', plotN, cpuTime, '-x');
xlabel('Matrix Size (n/100)');
ylabel('Cpu Time, s');
legend('Cn^3','CpuTime');