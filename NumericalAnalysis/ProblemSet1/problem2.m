x = [10e0 10e-1 10e-2 10e-3 10e-4 10e-5 10e-6 10e-7 10e-8 10e-9];
y = sin(x.*x)./(x.*x);
z = (sin(x).*sin(x))./(x.*x);

y'
z'