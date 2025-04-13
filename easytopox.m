
classdef easytopox <hgsetget
    properties
        BrainField
        cmap
        BrainSurface
        map
        img
        mni
        df
        colors
        NV
        NF
        neiber
        label
        facearea
        v2face
        pcutoff
        cutoff
        pn
        clustersize
        t
        p
        n
        minN
        minP
        std
        mean
    end % end of properties
    methods
        function obj=setParameter(obj,pcutoff,df,pn,clustersize,minN,minP)
            % pcutoff [0.01 0.001];
            %df Nsubj-1;
            %pn render [-1 1];
            %minN if a point n data < ignor.
            %minP if a ROI prob <  ignor in report
            obj.pcutoff=pcutoff;
            obj.cutoff=abs(tinv(pcutoff,df)); % a crude p level.  in ttest p for each point will e provided
            obj.pn=pn;
            obj.df=df;
            obj.clustersize=clustersize;
            obj.minN=minN;
            obj.minP=minP;
            obj.p=[];
        end
        function clusterStat=statCluster(obj)
            clusterStat=[];
            for pn=1:length(obj.pn)
                ml=max(obj.label{pn});
                
                vt=obj.label{pn};
                nn=0;
                clusterStat(pn).label=vt;
                for i=1:ml
                    ind=find(vt==i);
                    if length(ind)<obj.clustersize; continue; end
                    nn=nn+1;
                    temp=obj.t(ind);
                    [mv,ii]=max(temp*obj.pn(pn));
                    clusterStat(pn).N(nn)=length(ind);
                    clusterStat(pn).peakt(nn)=mv*obj.pn(pn);
                    j=ind(ii);
                    try
                    clusterStat(pn).peakxyz(nn,:)=obj.BrainSurface.vertices(j,:);
                    clusterStat(pn).ind(nn)=j;
                    clusterStat(pn).peakp(nn)=obj.p(j);
                    end
                    try
                    clusterStat(pn).peakbeta(nn)=obj.mean(j);
                    if obj.mean(j)*mv<0;
                        aaa=0
                    end
                    save('ojb','obj');
                    clusterStat(pn).obj=load('ojb');;
                    clusterStat(pn).peakerror(nn)=obj.std(j)/sqrt(obj.n(j));
                    clusterStat(pn).n(nn)=obj.n(j);
                    end
                    facel=[];
                    for ii=1:length(ind)
                        iii=ind(ii);
                        facel=[facel obj.v2face{iii}];
                    end
                    clusterStat(pn).facel{nn}=facel;
                    clusterStat(pn).area(nn)=sum(obj.facearea(facel));
                end
                try
                    clusterStat(pn).result7=nfri_anatomlabel_modifiedX(clusterStat(pn).peakxyz, '', 10, 7);
                    clusterStat(pn).result4=nfri_anatomlabel_modifiedX(clusterStat(pn).peakxyz, '', 10, 4);
                    clusterStat(pn).result5=nfri_anatomlabel_modifiedX(clusterStat(pn).peakxyz, '', 10, 5);
                end
            end
        end
        function obj=ttest(obj,mall)
            N=37694;
            obj.t=zeros(N,1)+NaN;
            obj.n=zeros(N,1);
            obj.p=ones(N,1);
            for i=1:N
                temp=mall(i,:);
                temp(temp==0)=[];
                temp(isnan(temp))=[];
                obj.n(i)=length(temp);
                if length(temp)<3;obj.t(i)=NaN;obj.p(i)=1; continue;end
                temp(isoutlier(temp))=[];
                [H,obj.p(i),CI,STATS] = ttest(temp) ;
                obj.t(i)=STATS.tstat;
                if obj.n(i)<obj.minN; obj.t(i)=0;obj.p(i)=1;continue;end
                obj.std(i)=STATS.sd;
                obj.mean(i)=mean(temp);
                if mod(i,300)==0; fprintf('.');end
            end
            %ind=find(p>0.5);
            %p(ind)=1-p(ind);
            try
            obj.p=obj.p/2;% so that make it one tail as spm
            obj.p(obj.n<obj.minN)=1;
            obj.t(obj.n<obj.minN)=0;
            fprintf('\n')
            catch
                'hi'
            end
        end
        function stat=ttestmap(obj,mall)
            %mostly useless. not to use
            for i=1:181
                for j=1:201
                    temp=mall(i,j,:);
                    temp(temp==0)=[];
                    temp(isnan(temp))=[];
                    if length(temp)==0; continue;end
                    [H,stat.p(i,j),CI,STATS] = ttest(temp) ;
                    stat.t(i,j)=STATS.tstat;
                    stat.n(i,j)=length(temp);
                end
                fprintf('.');
            end
            fprintf('\n');
        end
        function regression(obj,mall,X,C)
            %mall is 37694 * nsample matrix
            %X is covariance nsample * n regressor remember the last regressor has to be 1
            %C is  contract matrix  n final contrast * n regressor ,    eye(nregressor) will be ok.
            % output t  first 1 is ttest, 2nd t is first contrast
            % first result is Main effect, second... is the first ...
            % covaiance
            N=37694;
            N=size(mall,1);
            nC=size(C,1);
            obj.t=zeros(N,nC+1)+NaN;
            obj.n=zeros(N,1);
            obj.p=zeros(N,nC+1)+1;
            for i=1:N
                XX=X;
                temp=mall(i,:);
                XX(isnan(temp)| temp==0,:)=[];
                temp(isnan(temp)| temp==0 )=[];
                obj.n(i)=length(temp);
                if obj.n(i)<obj.minN; obj.p(i,:)=1; obj.t(i,:)=0;continue;end
                %if length(temp)==0;obj.t(i)=NaN; continue;end
                [H,obj.p(i,1),CI,STATS] = ttest(temp) ;
                XX(:,end+1)=1;
                Result=nirsRegX(temp',XX,C);
                obj.t(i,1)=STATS.tstat;
                obj.t(i,2:end)=Result.t(1,:);
                obj.p(i,2:end)=Result.p(1,:);
                if mod(i,300)==0; fprintf('.');end
            end
        end
        function LoadV(a,t,vol)
            a.img.val=zeros(37694,1)+NaN;
            pos = find(abs(t)>0);
            clear xyzcor
            [xyzcor(:,1) xyzcor(:,2) xyzcor(:,3)] = ind2sub(size(t), pos);
            xyz=cor2mni(xyzcor,vol.mat);
            ds=pdist2(a.BrainSurface.vertices,xyz);
            [d ind]=min(ds');
            mz=min(xyzcor(:,3));
            for i=1:length(ind);
                j=ind(i);
                if d(i)<6
                    a.img.val(i)=t(xyzcor(j,1),xyzcor(j,2),xyzcor(j,3));
                end
            end
        end
        function showmask(a,fn,cl,r) % fn is {'a' 'b'} cl  is [ -1 1 2]  -1 blue 1 2 some red yellow
            % r is the inflation radius
            for m=1:length(fn)
                vol=spm_vol(sprintf('%s.img',fn{m}));
                mask(:,:,:,m)=spm_read_vols(vol);
            end
            a.img.val=zeros(37694,1);
            for m=1:length(fn)
                pos = find(mask(:,:,:,m)>=0.5);
                clear xyzcor
                [xyzcor(:,1) xyzcor(:,2) xyzcor(:,3)] = ind2sub(size(mask(:,:,:,m)), pos);
                xyz=cor2mni(xyzcor,vol.mat);
                ds=pdist2(a.BrainSurface.vertices,xyz);
                ind=min(ds')<5;
                a.img.val(ind)=cl(m);
            end
            figure;
            a.setParameter([0.5 0.01],1,[1 -1],1,1,1)
            a.cutoff=[0.5 max(abs(cl)+0.5)];
            a.show3D();
            view(-90,1);
        end
        function  [h,hl] =showCluster(obj,clustersize,pn0)
            val=obj.img.val;
            obj.img.val=val*0;
            mm=-1;
            for pn=1:length(pn0)
                t0=pn0(pn);
                ii=find(obj.pn==t0);
                if isempty(ii);continue;end
                ml=max(obj.label{ii});
                vt=obj.label{ii};
                for i=1:ml
                    ind=find(vt==i);
                    Nl(i)=length(ind);
                    if Nl(i)<obj.clustersize; continue; end
                    obj.img.val(ind)=val(ind);
                end
            end
            [h,hl] =obj.show3D();
            obj.img.val=val;
        end
        
        function showClusterLabel(obj,clustersize,pn0)
            val=obj.img.val;
            obj.img.val=val*0;
            mm=-1;
            for pn=1:length(pn0)
                ml=max(obj.label{pn});
                vt=obj.label{pn};
                for i=1:ml
                    ind=find(obj.label{pn}==i);
                    Nl(i)=length(ind);
                    if Nl(i)<clustersize; vt(vt==i)=0; end
                end
                obj.img.val=obj.img.val+pn0(pn)*vt;
                if ml>mm; mm=ml;end
            end
            obj.show3D();
        end
        function cluster(obj)
            global neiborList val0 label1 l
            neiborList=obj.neiber;
            obj.p(isnan(obj.p))=1;
            for ii=1:length(obj.pn)
                pn=obj.pn(ii);
                if pn==0; continue;end
                v=obj.img.val*pn;
                if ~ isempty(obj.p)
                    p=obj.p;
                    p(v<=0)=1;
                    val0=p<obj.pcutoff(1);  % if you just did obj.ttest
                else
                    val0=v>obj.cutoff(1);
                end
                label1=val0*0;
                l=0;
                for i=1:obj.NV
                    if val0(i)>0; l=l+1; SearchCluster(i); end
                    if mod(i,100)==0; fprintf('.');end
                end
                fprintf('\n');
                obj.label{ii}=label1;
            end
        end
        function obj=preparecluster(obj)
            obj.NV=size(obj.BrainSurface.vertices,1);
            obj.NF=size(obj.BrainSurface.faces,1);
            for i=1:obj.NV
                obj.neiber{i}=[];
                obj.v2face{i}=[];
            end
            for j=1:obj.NF
                temp=obj.BrainSurface.faces(j,:);
                for ii=1:3
                    temp1=temp;
                    i=temp(ii);
                    temp1(ii)=[];
                    obj.neiber{i}=[obj.neiber{i} temp1];
                    obj.v2face{i}=[obj.v2face{i} j];
                end
                obj.facearea(j)=triangleP(obj.BrainSurface.vertices(temp,:));
                if mod(j,1000)==0; fprintf('.');end
            end
            fprintf('\n');
            for i=1:obj.NV
                obj.neiber{i}=unique( obj.neiber{i});
                obj.v2face{i}=unique( obj.v2face{i});
            end
        end
        function obj = easytopox()
            load('BrainField.mat');
            obj.BrainField = BrainField;
            clear BrainField;
            obj.BrainSurface = load('BrainSurface.mat');
            load cmap_jet;
            obj.cmap = cmap;
            load('rendercolors.mat');
            obj.colors=colors;
        end
        function obj = setmni(obj,mni)
            mni(:,4)=1+(mni(:,1)>0);
            obj.mni = mni;
            obj.loadT(mni(:,1)); % set up phi thea
        end
        function obj = loadT(obj,val)
            try
            [theta, phi, radius, value, region] = topo_map_rot(obj.mni, val);
            obj.map.val = value;
            obj.map.theta = theta;
            obj.map.phi = phi;
            end
            obj.img.val = topo_surf_rot(obj.mni, val, obj.BrainSurface.vertices);
            
        end
        function h = showphitheta(obj)
            map = obj.map.val;
            if  ~isempty(obj.cutoff) & abs(diff(obj.cutoff))>0
                h = imagesc(obj.map.phi(1,:),obj.map.theta(:,1),map,[-cutoff1 cutoff1]); set(gca,'Ydir','normal');
            else
                h = imagesc(obj.map.phi(1,:),obj.map.theta(:,1),map,[-max(abs(map(:))) max(abs(map(:)))]); set(gca,'Ydir','normal');
            end
            hold on; contour(obj.map.phi,obj.map.theta,obj.BrainField,5,'k-','linewidth',1);
            xlabel('\phi (degree)'); ylabel('\theta (degree)');
            colormap(obj.cmap); colorbar;
        end
        function [h,hl] = show3D(obj)
            vertvalue = obj.img.val;
            vertcolor=val2color(vertvalue,obj.cutoff,obj.pn,obj.colors);
            [h,hl] = plotimageeasy(obj.BrainSurface.vertices,obj.BrainSurface.faces, vertcolor, [0 0 1]);
            colormap(obj.cmap);
        end
        function h=plotchannel(obj,plotBrain,plotText)
            h = plotoptodes([],obj.BrainSurface.faces,obj.mni,4,[0 0 1]);
            %h = plotoptodes(obj.BrainSurface.vertices,obj.BrainSurface.faces,obj.mni,3,[0 0 1]);
            if plotBrain>0
                val=obj.img.val;
                obj.img.val=0*obj.img.val;
                obj.show3D();
                obj.img.val=val;
            end
            if plotText==0; return;end
            d=1.1;
            for i=1:size(obj.mni,1);
                text(obj.mni(i,1)*d,obj.mni(i,2)*d,obj.mni(i,3)*d,int2str(i));
            end
        end
        function reportCluster(obj,c,str,fid);
            fprintf(fid,',,,,,,,,,%s\n',str);
            fprintf(fid,'T,');
            fprintf(fid,'p,');
            %fprintf(fid,'N pixel %d,',N(i));
            fprintf(fid,'X,Y,Z,');
            fprintf(fid,'df,');
            fprintf(fid,'CM2\n');

            for pn=1:length(obj.pn);
                cc=c(pn);
                %fprintf(fid,'part %d\n',pn);
                N=cc.N(cc.N>=obj.clustersize);
                try
                    for i=1:length(cc.peakt);
                        df=cc.n(i)-1;
                        if df<obj.minN-1; continue;end
                        if cc.peakp(i)>obj.pcutoff(1); continue;end
                        fprintf(fid,'%.2f,',cc.peakt(i));
                        fprintf(fid,'%.2f,',cc.peakp(i)*100);
                        %fprintf(fid,'N pixel %d,',N(i));
                        fprintf(fid,'%.1f, %.1f,%.1f,',cc.peakxyz(i,1),cc.peakxyz(i,2),cc.peakxyz(i,3));
                        fprintf(fid,'%d,',cc.n(i)-1);
                        fprintf(fid,'%.1f,',cc.area(i)/100);
                        for jj = [4,5,7]
                          eval(sprintf('r=cc.result%d{i};',jj))
                        
                        for k=1:length(r.Perc);
                            p=r.Perc(k);
                            if p<obj.minP; continue;end
                            if k==1
                                fprintf(fid,'%s,%.2f\n',strrep(r.ROIlist{k},',',';'),p);
                            else
                                fprintf(fid,',,,,,,,%s,%.2f\n',strrep(r.ROIlist{k},',',';'),p);
                            end
                        end
                        end
                    end
                catch
                    a=1;
                end
            end
        end
         function rendereasytopo4(obj,titlestr,fn,pn0,plotseparate);
             % assum obj.img.val is filled
            if nargin==4
                plotseparate=1;
            end
            i=1;
            lightv(i).Color=0.8*[ 1 1 1];
            lightv(i).Style='local';
            lightv(i).Position=[0.9931 -15.9641 1.3828e+03];
            lightv(i).Visible='on';
            i=2;
            lightv(i).Color=0.8*[ 1 1 1];
            lightv(i).Style='infinite';
            lightv(i).Position=[0 -1 1];
            lightv(i).Visible='on';
            i=3;
            lightv(i).Color=0.8*[ 1 1 1];
            lightv(i).Style='infinite';
            lightv(i).Position=[0.7500 0 0.5000];
            lightv(i).Visible='on';
            i=4;
            lightv(i).Color=0.8*[ 1 1 1];
            lightv(i).Style='infinite';
            lightv(i).Position=[-0.7500 0 0.5000];
            lightv(i).Visible='on';
            if length(pn0)==2; 
                extn='PN';
            else
                if pn0(1)>0
                    extn='P';
                else
                    extn='N';
                end
            end
            extn=sprintf('%s%d',extn,round(obj.pcutoff(1)*100000));
            plotseparate=0;
            figure;%h2=figure;
            vl={[0 90] [0 45] [90 0] [270 0]};
            ll={'headlight' 'right' 'right' 'left'};
            oname={'top' 'right' 'back' 'right' };
            [a,b]=fileparts(fn);
            b=strrep(b,',','');
            fd=fullfile(strrep(a,',',''),'img_easytopo');
            if ~exist(fd,'dir'); mkdir(fd);end
            if plotseparate==0;
                ha=tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);
            end
            for vv=4:-1:1
                if plotseparate==0;
                    axes(ha(vv));
                end
                [~,h] =obj.show3D;
                view(vl{vv}(1),vl{vv}(2));
                %camlight(h,ll{vv});
                set(h,'Position',lightv(vv).Position,'Color',lightv(vv).Color,'Style',lightv(vv).Style);
                if plotseparate==1
                    saveas(gca,sprintf('%s/%s_%s_%s.png',fd,b,oname{vv},extn));
                end
            end
            if plotseparate==1 ;return;end
            saveas(gcf,sprintf('%s/%s_%s.png',fd,b,extn));
        end
        function rendereasytopo(obj,titlestr,fn,clustersize,pn0,plotseparate);
            if nargin==5
                plotseparate=1;
            end
            
            FourPanel=1;
            LRP=2;
            RP=3;
            LP=4;
            AP=5;

            i=1;
            lightv(i).Color=0.8*[ 1 1 1];
            lightv(i).Style='local';
            lightv(i).Position=[0.9931 -15.9641 1.3828e+03];
            lightv(i).Visible='on';
            i=2;
            lightv(i).Color=0.8*[ 1 1 1];
            lightv(i).Style='infinite';
            lightv(i).Position=[0 -1 1];
            lightv(i).Visible='on';
            i=3;
            lightv(i).Color=0.8*[ 1 1 1];
            lightv(i).Style='infinite';
            lightv(i).Position=[0.7500 0 0.5000];
            lightv(i).Visible='on';
            i=4;
            lightv(i).Color=0.8*[ 1 1 1];
            lightv(i).Style='infinite';
            lightv(i).Position=[-0.7500 0 0.5000];
            lightv(i).Visible='on';
            if length(pn0)==2; 
                extn='PN';
            else
                if pn0(1)>0
                    extn='P';
                else
                    extn='N';
                end
            end
            extn=sprintf('%s%d',extn,round(obj.pcutoff(1)*100000));
            plotseparate=0;
            figure;%h2=figure;
            vl={[0 90] [0 45] [90 0] [270 0]};
            ll={'headlight' 'right' 'right' 'left'};
            oname={'top' 'right' 'back' 'right' };
            [a,b]=fileparts(fn);
            b=strrep(b,',','');
            fd=fullfile(strrep(a,',',''),'img_easytopo');
            if ~exist(fd,'dir'); mkdir(fd);end
            if plotseparate==0;
              ha=tight_subplot(2,2,[.01 .03],[.1 .01],[.01 .01]);
            end
            for vv=4:-1:1
                if plotseparate==0;
                    axes(ha(vv));
                end
                [~,h] =obj.showCluster(clustersize,pn0);
                view(vl{vv}(1),vl{vv}(2));
                %camlight(h,ll{vv});
                set(h,'Position',lightv(vv).Position,'Color',lightv(vv).Color,'Style',lightv(vv).Style);
                if plotseparate==1
                    saveas(gca,sprintf('%s/%s_%s_%s.png',fd,b,oname{vv},extn));
                end
            end
            if plotseparate==1 ;return;end
            saveas(gcf,sprintf('%s/%s_%s.png',fd,b,extn));
        end
    end
end
