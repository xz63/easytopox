function facecolor=val2color(v1,cutoff,pn0,colors);
threshold=cutoff(1);
topZ=cutoff(2);
facecolor=ones(length(v1),3);
for pn=-1:2:1
    if isempty(find(pn==pn0));continue;end
    intensity = pn*v1;
    if pn==-1
	mymap=colors{1};
	mymap=mymap(end:-1:1,:);
    else
	mymap=colors{2};
    end
    ii=find(intensity>threshold);
    intensity1=intensity(ii);
    intensity1(intensity1>topZ)=topZ;
    cindex=floor(99*(intensity1-threshold)/(topZ-threshold))+1;
    facecolor(ii,:)=mymap(cindex,:);
end
return
ind=1:10;
for i=1:3
im(:,:,i)=reshape(facecolor(:,i),sz);
subplot(2,2,i+1);
imagesc(im(ind,ind,i));
end
figure;
cl='rgb'
for i=1:3
subplot(2,2,i);
plot(val(:),facecolor(:,i),cl(i));
end

