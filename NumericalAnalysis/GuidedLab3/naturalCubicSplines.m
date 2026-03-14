function [A,rvec,c,a,b,d] = naturalCubicSplines(x,y)
    if length(x) < 2 | length(y) < 2
        error('Input vectors are not long enough, please make them of length two or more and be of the same size')
    end
    if length(x) ~= length(y)
        error('Input vectors are not the same length')
    end
    A=zeros(length(x));
    rvec=zeros(length(x),1);
    h=zeros(1,length(x));

    % Begin Problem 1
    A(1,1) = 1;
    A(length(x),length(x)) = 1;

    h(1) = x(2) - x(1);

    for i = 2:length(x)-1
        h(i) = x(i+1) - x(i);
        A(i,i-1) = h(i-1);
        A(i,i) = 2*(h(i-1)+h(i));
        A(i,i+1) = h(i);
        rvec(i) = 3*((y(i+1)-y(i))/h(i) - (y(i)-y(i-1))/h(i-1));
    end
    % End Problem 1

    % Problem 2 - Solve the system of equations
    c = A \ rvec;

    % Problem 3 - Calculate the coefficients a, b, and d
    a = zeros(length(x),1);
    for i=1:length(x)
        a(i) = y(i);
    end

    b = zeros(length(x)-1, 1);
    d = zeros(length(x)-1, 1);

    for i = 1:length(x)-1
        b(i) = (c(i+1) - c(i))/(3*h(i));
        d(i) = (a(i+1) - a(i))/h(i) - (h(i)*(2*c(i) + c(i+1)))/3;
    end

    % Problem 4 - use the given code to perform
    coeffs=[];
    for k=1:length(x)-1
        coeffs(k,:)=[d(k) c(k) b(k) a(k)];
    end

    breaks=[x]
    pp=mkpp(breaks,coeffs)
    ystar=ppval(pp,xstar)

end