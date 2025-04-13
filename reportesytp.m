%function h=report(obj,clusterStat,cname,fd,fn)
if 1>0
    clustersize=20;
    pn=[-1 1];
    df=120;
    pcutoff=[0.05 0.005];
    obj.setParameter(pcutoff,df,pn,clustersize);
    obj=preparecluster();
    fd='result';
    cname='test';
    fn='a.html';
    if ~exist(fd,'dir'); mkdir(fd);end
    h=fn;
    s=class(fn);
    if s(1)=='c'; h=htmlTable(fullfile(fd,fn)); end
    fn=fullfile(fd,cname);
    figure;
    [~,hc]=obj.show3D();
    view(0,90); camlight(hc,'headlight'); refresh; saveas(gca,[fn 'top.png'],'png');
    view(90,0); camlight(hc,'right'); refresh; saveas(gca,[fn 'right.png'],'png');
    view(270,0); camlight(hc,'left'); refresh; saveas(gca,[fn 'left.png'],'png');
    view(0,0); camlight(hc,'headlight'); refresh; saveas(gca,[fn 'back.png'],'png');
    view(180,15); camlight(hc,'headlight'); refresh; saveas(gca,[fn 'front.png'],'png');
    rnds={'right' 'left' 'top' 'front' 'back' ''};
    cn=cname;
    pValue =obj.pcutoff(1);
    obj.cluster();
    stat=obj.ttest(allbetaimg(:,:,1));
    clusterStat=obj.statCluster(stat);
end % done
if length(obj.pn)>1;
    EXT='';
else
    if obj.pn(1)==-1
        EXT='Negative';
    else
        EXT='Positive';
    end
end
ii=1;
clear tt
tt{ii}=sprintf('@@%s,  Cluster size threshold: %d P threshold %f' ,datestr(now,'dd-mmmm-yy'), obj.clustersize,obj.pcutoff(1)); ii=ii+1;
tt{ii}=sprintf('<p style="page-break-before: always">\n GLM Contrast: %s ',cname); ii=ii+1;
clear imgtable
for jj=1:3
    fnim=sprintf('%s%s.png',cn,rnds{jj});
    imgtable{1,jj}=fnim;
end
h.addImgTable(imgtable,'',440,326,rnds(1:3));
clear imgtable
for jj=1:2
    fnim=sprintf('%s%s.png',cn,rnds{3+jj});
    imgtable{1,jj}=fnim;
end
imgtable{1,3}='';
h.addImgTable(imgtable,'',400,300,rnds(3+[1:3]));
i=1;
for pn=1:length(obj.pn)  
    if isempty(clusterStat(pn).N); imgtable{1,i}='No supra-threshold cluster';i=i+1;continue;end
    for j=1:length(clusterStat(pn).N)
       fprintf(obj.fid,'Peak T: %.2f, P: %.5f , DF: %d, N Voxel %d',clusterStat(pn).peakt(j), ...
          clusterStat(pn).peakp(j), clusterStat(pn).n(j));
        fprintf(obj.fid,'Peak MNI coordinate %f,',clusterStat(pn).peakxyz(j,:));
        temp=clusterStat(pn).result7{j};
        for ii=1:length(clusterStat(pn).result7{j}.Perc)
            if temp{1}.Perc(ii)<.1; break; end
            fprintf('%s %f\n',temp{1}.ROIlist{ii},temp{1}.Perc(ii) );
            res{mm}.result{jj}=  sprintf('%s %f\n',temp{1}.ROIlist{ii},temp{1}.Perc(ii) ); jj=jj+1;
        end
    end
       
    res={};
    result={};
    if ~exist(fn);result{1}='No supra-threshold cluster'; continue; end
    fprintf('%s exist\n', fn); load(fn);
    fn=sprintf('/%s/%s%s',x.procesName{R+1},cn,x.signalType{k});
    
end
n=0;
n=n+1;
result{n}=' ';
for jj=1:length(res{i}.result)
    n=n+1; result{n}=res{i}.result{jj};
end

h.close();
