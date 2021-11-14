function [ HFEN_m ] = hfen( U,ref )
%HFEN High-Frequency Error

U=abs(U); r=abs(ref);

for t = 1:size(r,3)
    alpha = sum(dot(U(:,:,t),r(:,:,t)))/(sum(dot(U(:,:,t),U(:,:,t))));
    U(:,:,t)=(alpha)*U(:,:,t);
    HFEN_m(t)=norm(imfilter(abs(U(:,:,t)),fspecial('log',15,1.5)) - imfilter(abs(r(:,:,t)),fspecial('log',15,1.5)),'fro')./norm( imfilter(abs(r(:,:,t)),fspecial('log',15,1.5)),'fro');
end

HFEN_m=mean(HFEN_m);

