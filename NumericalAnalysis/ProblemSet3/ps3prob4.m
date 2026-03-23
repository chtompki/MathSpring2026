xrng=[-1,3]; %interval
n=5; %n intervals, n+1 data points

funInterp=@(x) x.*exp(-x); %define function to approximate
funInterpPrime=@(x) exp(-x)+x.*exp(-x); %define derivative

%% polynomial interpolations
%uniformly spaced points
xunif=linspace(xrng(1),xrng(2),n+1);
yunif=funInterp(xunif);
%chebyshev points
xcheb=1+2*cos((2*(0:n)+1)/(2*(n))*pi );
ycheb=funInterp(xcheb);
%legendre points
syms xs %use matlab symbolic toolbox
xlgndrS=vpasolve(legendreP(n+1,xs) == 0); %symbolically solve for roots
xlgndr=double(xlgndrS); %convert to double precision
ylgndr=funInterp(xlgndr);


xP=linspace(xrng(1),xrng(2),101); %x values for plotting
yP=funInterp(xP); %actual function at ypoints

%use polyfit with degree n to get nth degree polyomial interpolant
acoeff=polyfit(xunif,yunif,n);%get coefficients
yPunif=polyval(acoeff,xP); %evaluate at plot points

acoeff=polyfit(xcheb,ycheb,n);%get coefficients
yPcheb=polyval(acoeff,xP);%evaluate at plot points

acoeff=polyfit(xlgndr,ylgndr,n);%get coefficients
yPlgndr=polyval(acoeff,xP);%evaluate at plot points


%plot interpolants
figure;hold on;
plot(xP,yP,'-')
plot(xP,yPunif,':')
plot(xP,yPcheb,'--')
plot(xP,yPlgndr,'-.')
legend('Standard Interpolation', 'Chebyshev', 'Legandre')

%plot interpolation error
figure;hold on;
plot(xP,yP-yPunif,':')
plot(xP,yP-yPcheb,'--')
plot(xP,yP-yPlgndr,'-.')
legend('Standard Interpolation', 'Chebyshev', 'Legandre')