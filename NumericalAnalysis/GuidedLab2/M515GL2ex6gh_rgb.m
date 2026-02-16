%%
rgbGH=imread('GraceHopper.jpg');
figure();
image(rgbGH), axis image;  %plot color image
[N,M,~]=size(rgbGH);

rGH=rgbGH(:,:,1);
gGH=rgbGH(:,:,2);
bGH=rgbGH(:,:,3);

bwGH=rgb2gray(rgbGH); %convert to grayscale
imGH=double(bwGH); %convert from unsigned integers to double for calculutions


%plt grayscale image
figure();
colormap(gray(256));
image(imGH);
daspect([1 1 1]) %this preserves aspect ratio, could also use "axis image"
title('Original');


%Compute singular value decomposition
[U S V]=svd(imGH);

[Ur Sr Vr]=svd(double(rGH));
[Ug Sg Vg]=svd(double(gGH));
[Ub Sb Vb]=svd(double(bGH));

%plot singular values on semilog scale.
%Notice how quickly the magnitude drops.
figure;
semilogy(diag(S))
ylabel('Singular Values')
xlabel('n')


%pick number of singular values to use for reconstructing image
Nsvals=[200, 100, 50, 25, 10, 5, 2, 1];

%plot "compressed" images by only include the largest ns singular values
for jj=1:length(Nsvals)
    ns=Nsvals(jj);
    
    imNs=U(:,1:ns)*S(1:ns,1:ns)*V(:,1:ns)';

    imNsRGB=zeros(N,M,3);
    imNsRGB(:,:,1)=Ur(:,1:ns)*Sr(1:ns,1:ns)*Vr(:,1:ns)';
    imNsRGB(:,:,2)=Ug(:,1:ns)*Sg(1:ns,1:ns)*Vg(:,1:ns)';
    imNsRGB(:,:,3)=Ub(:,1:ns)*Sb(1:ns,1:ns)*Vb(:,1:ns)';


    figure();
    colormap(gray(256));
    image(imNs);
    daspect([1 1 1])
    title([num2str(ns)  ' singular values']);

    figure();
    image(uint8(imNsRGB));
    daspect([1 1 1])
    title([num2str(ns)  ' singular values']);
    
end