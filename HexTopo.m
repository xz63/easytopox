classdef HexTopo <easytopox
    properties
        distM
        idx
        center
        surfacehex256
        ii
    end % end of properties
    methods
        function obj=HexTopo()
            obj=obj@easytopox();
            load surfacehex256 %surfacehex256 center256
            obj.idx=idx;
            obj.center=center256;
            obj.surfacehex256=surfacehex256;
            distM=squareform(pdist(center256));
            [~,obj.ii]=sort(distM)
        end
        function obj=mkHex(obj)
            load mask
            obj.LoadV(mask,vol);
            ind=find(~isnan(obj.img.val));
            nC=256;
            [idx,center256]=kmeans(obj.BrainSurface.vertices(ind,:),nC);
            obj.img.val(ind)=3*(idx-(nC/2))/(nC/2);
            obj.show3D();
            obj.img.val(ind)=idx;
            surfacehex256=obj.img.val;
            save surfacehex256 surfacehex256 center256 idx
            obj.idx=idx;
            obj.center=center256;
            obj.surfacehex256=surfacehex256;
            distM=squareform(pdist(center256));
            [~,ii]=sort(distM)
        end
        function obj=demoHex(obj)
            subidx=obj.surfacehex256;
            for i=1:7 %pick 7 neighor,
                subidx(subidx==obj.ii(i,7))=3000; % the 7th patch as center
            end
            subidx=subidx==3000;
            obj.img.val=3*subidx;
            obj.show3D();
        end
        function corimg=spatial_correlation(obj,val) % all vertex by 2
            for n=1:256
                subidx=obj.surfacehex256;
                for i=1:7 %pick 7 neighor,
                    subidx(subidx==obj.ii(i,n))=3000; % the 7th patch as center
                end
                temp=val(subidx==3000,:,:);
                tt=sum(temp(:,:),2);
                temp(isnan(tt),:,:)=[];
                if isempty(temp)
                    corv(n)=NaN;
                else
		    t1=temp(:,:,1);
		    t2=temp(:,:,2);
                    tt=corrcoef(t1(:),t2(:));
                    corv(n)=tt(2,1);
                end
            end
            obj.setmni(obj.center);
            obj.loadT(corv);
            corimg=obj.img.val;
        end
    end
end

