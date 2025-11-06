%Calculate the diagonal elements of F'F externally
%Obtain the F'F diagonal element matrices under different ker_std
%Since the matrix is only related to the image size and the two sample images are of the same size, it is acceptable to use either sample image
clc,clear

I=imread('cameraman_big.jpg');I = rgb2gray(I);X=im2double(I);X=imresize(X,0.5,'nearest');
%I = imread('house.png');X=im2double(I);
[m,n]=size(X);
ker_size = 15; 
%%
ker_std =0.5;
p = fspecial('gaussian', ker_size, ker_std);
K = imfilter(X,p,'symmetric');
F_tF=zeros(m,n);
    for i=1:m
        for j=1:n
            E=zeros(m,n);
            E(i,j)=1;
            Temp = idct2(E);
            Temp1 = imfilter(Temp,p,"symmetric");
            F_tF(i,j)=norm(Temp1,'fro')^2;
        end
    end
save('F_tF_small','F_tF')
%%
ker_std =1;
p = fspecial('gaussian', ker_size, ker_std);
K = imfilter(X,p,'symmetric');
F_tF=zeros(m,n);
    for i=1:m
        for j=1:n
            E=zeros(m,n);
            E(i,j)=1;
            Temp = idct2(E);
            Temp1 = imfilter(Temp,p,"symmetric");
            F_tF(i,j)=norm(Temp1,'fro')^2;
        end
    end
save('F_tF_mid','F_tF')
%%
ker_std =1.5;
p = fspecial('gaussian', ker_size, ker_std);
K = imfilter(X,p,'symmetric');
F_tF=zeros(m,n);
    for i=1:m
        for j=1:n
            E=zeros(m,n);
            E(i,j)=1;
            Temp = idct2(E);
            Temp1 = imfilter(Temp,p,"symmetric");
            F_tF(i,j)=norm(Temp1,'fro')^2;
        end
    end
save('F_tF_big','F_tF')