time=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24];
temp=[58 58 58 58 57 57 57 58 60 64 67 68 66 66 65 64 63 63 62 61 60 60 59 58];

xstar = zeros(100);
for i = 1:100
    xstar(i) = i*0.24;
end

ystar = naturalCubicSplines(time,temp,xstar);

figure;
hold on;
plot(xstar,ystar, 'k-');
plot(time,temp, 'ko');
xlabel('Time (hours)');
ylabel('Temperature (°F)');
hold off;