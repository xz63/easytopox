function reportCluster(obj,c,str,fid,minn,minp);
%minn=4; minp=0.25;
%fid=fopen([fn '.txt'],'w');
%c= clusterStat{2,2};
fprintf(fid,'%s\n',str);
for pn=1:length(c);
    cc=c(pn);
    fprintf(fid,'part %d\n',pn);
    N=cc.N(cc.N>=obj.clustersize);
    for i=1:length(cc.peakt);
	  df=cc.n(i)-1;
	  fprintf(fid,'N pixel %d,',N(i));
	  fprintf(fid,'MNI %.1f, %.1f,%.1f,',cc.peakxyz(i,1),cc.peakxyz(i,2),cc.peakxyz(i,3));
	  fprintf(fid,'T %.2f,',cc.peakt(i));
	  fprintf(fid,'p %.2f,',cc.peakp(i)*100); 
	  fprintf(fid,'df %d,',cc.n(i)-1); 
	  fprintf(fid,'area %.2f mm\n',cc.area(i)); 
	  r=cc.result7;
	  for j=1:length(r)
	      for k=1:length(r{j}.Perc);
		  p=r{j}.Perc(k);
		  if p<minp; continue;end
		  fprintf(fid,'\t%s,%.2f\n',r{j}.ROIlist{k},p);
	      end
	  end
      end
  end
