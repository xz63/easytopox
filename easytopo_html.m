classdef easytopo_html <easytopox
    properties
        h
	fid
    end % end of properties
    methods
        function fid=starthtml(obj,fn)
	    obj.h=htmlTable(fullfile(htmlfd,fn));
	    fid=obj.h.fid;
	    fprintf(fid,'Threshold P %f, cluster threshold %d, date: %s\n', obj.pcutoff(1), obj.clustersize,datestr(now));
	    fprintf(fid,'<p style="page-break-before: always">\n'); 
	end
	function closehtml(obj);
	    obj.h.close();
	end
	function addimage(obj,titlestr,fn)
	    ll={'headlight' 'right' 'left'};
	    fprintf(fid,'%s <br>',titlestr);
	    fprintf(fid,'<table border="0">\n');
	    fprintf(fid,'<tr>\n');
	    w2=571;
	    h2=193;
	    for c=1:3
		fprintf(fid,'<img src="./imgs/%s_%s.png" width="%d" height="%d" >\n',fd,fn,ll{c},w2,h2);
	    end
	    fprintf(fid,'</tr>\n');
	    fprintf(fid,'</table>\n');
	end
	function addTable(obj)
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
                            if p<obj.minP; continue;end
                            fprintf(fid,'\t%s,%.2f\n',r{j}.ROIlist{k},p);
                        end
                    end
                end
            end
        end
        function rendereasytopo(obj,titlestr,fn,clustersize,pn0);
            figure;
            vl={[0 90] [90 0] [270 0]};
            ll={'headlight' 'right' 'left'};
	    if ~exist('img','dir'); mkdir('img');end
            for vv=3:-1:1
                subplot(2,2,vv);
                [~,h] =obj.showCluster(clustersize,pn0);
                view(vl{vv}(1),vl{vv}(2));
                camlight(h,ll{vv});
		saveas(gca,sprintf('img/%s_%s.png',fn,ll{vv}));
            end
            %title(strrep(titlestr,'_',' '));
            %saveas(gcf,sprintf('%s.pdf',fn));
        end
    end
end
	end
    end
end
