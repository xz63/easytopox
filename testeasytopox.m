obj=easytopox();
obj=obj.preparecluster;
load test.mat
for i=1:8
obj.setmni(xyz(:,:,i));
obj.loadT(beta(i,:)');
temp(i,:)=obj.img.val;
end
pcutoff=[0.05 0.005]; % the threshold and the upper bound;
df=8-1; % df
pn=[-1 1];   % both blue the red  [1] positive only [-1] negative only
clustersize=20; % cluster threhold
minN=5; %  not all points has data for every one, because head sizes are different
% at the there are minN subject, the T value is valid
minP  =0.25;  % for identify a point as a part of anatomicl location
obj.setParameter(pcutoff,df,pn,clustersize,minN,minP);
figure;
obj=obj.ttest(temp');
obj.img.val=obj.t;
obj.cluster();
 clusterStat=obj.statCluster();
obj.showCluster(obj.clustersize,obj.pn);

