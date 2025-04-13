pn0=[-1];
cutoff=[0.3 0.8];
load('rendercolors.mat');
ind=1:10;
val=zeros(11,11);
clear im
for i=1:11
    val(:,i)=-1:.2:1;
end
figure;
subplot(2,2,1); imagesc(val);
sz=size(val);
v1=val(:);
facecolor=val2color(v1,cutoff,pn0,colors);
for i=1:3
im(:,:,i)=reshape(facecolor(:,i),sz);
end
subplot(2,1,2);
imshow(im);
