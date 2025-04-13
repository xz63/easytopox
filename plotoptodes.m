function h = plotoptodes(vertices, faces, mni, radius, viewangle)
    h=-1;
    %plot brain surface
    if ~isempty(vertices)
    vertcolor = repmat([1 1 1], length(vertices), 1);
    h = plotsurf(vertices,faces,'facealpha',.3,'edgealpha',0,'facecolor','interp');
    set(h,'FaceVertexCData',vertcolor, 'SpecularStrength',0, 'AmbientStrength',0.4);  daspect([1 1 1]); axis off;
    set(h,'alpha',0.5)
    view(3); view(viewangle); 
    camlight headlight; lighting gouraud;
    clear vertcolor
    end

    hold on;
    %plot optodes
    pts = mni;
    rads = repmat(radius,size(pts,1),1);
    cols = repmat([1 0 0],size(pts,1),1);
    plotbubble3D(pts,rads,cols);


