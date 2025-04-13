function a=SearchCluster(i)
global val0 label1 l neiborList
val0(i)=0;
label1(i)=l;
temp=neiborList{i};
n=0;
for j=1:length(temp);
    k=temp(j);
    if val0(k)>0; n=SearchCluster(k)+n; end
end
a=n>0;
