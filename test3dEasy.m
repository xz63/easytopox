addpath /media/BFL/BFL-Share/BFLFNIR/ProgramTools/EasyTopoX
a=img;
for i=1:91;
    a(i,:,:)=3*(i-45)/91;
end
vol=spm_vol('/media/BFL/BFL-Share/BFLFNIR/Data/Haskins_treatment/spike_removeX/All/Clean/VisualR/mask.img'); 
mask=spm_read_vols(vol);
a=a.*mask;
obj=easytopox;
obj.setParameter([0.50 0.1],10,[-1 1],1)
obj.LoadV(a,vol);
obj.show3D();
ind=find(~isnan(obj.img.val));
nC=256;
[idx,c]=kmeans(obj.BrainSurface.vertices(ind,:),nC);
obj.img.val(ind)=3*(idx-(nC/2))/(nC/2);
obj.show3D();
obj.img.val(ind)=idx;
surfacehex256=obj.img.val;
center256=c;
save surfacehex256 surfacehex256 center256

d=squareform(pdist(center256));
[~,ii]=sort(d)[~,ii]=sort(d);
ii(1:10,8);
the 8th hexgon 1:10 nearest, first one is 8 useless.
subidx=surfacehex256;
for i=1:7 %pick 7 neighor,  
subidx(subidx==ii(i,7))=3000; % the 7th patch as center
end
subidx=subidx==3000;
obj.img.val=3*subidx;
obj.show3D();


