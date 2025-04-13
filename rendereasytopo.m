function rendereasytopo(obj,titlestr,fn);
figure;
vl={[0 90] [90 0] [270 0]};
ll={'headlight' 'right' 'left'};
%obj.ccluster();
for vv=3:-1:1
    subplot(2,2,vv);
    [~,h]= obj.show3D();
    view(vl{vv}(1),vl{vv}(2));
    camlight(h,ll{vv});
end
title(strrep(titlestr,'_',' '));
saveas(gcf,sprintf('%s.pdf',fn));
end

