function [h ,hl] = plotimage(vertices, faces, vertcolor, viewangle)

    h = plotsurf(vertices, faces, 'edgealpha', 0, 'facecolor', 'interp');
    set(h,'FaceVertexCData',vertcolor, 'SpecularStrength',0, 'AmbientStrength',0.4);  daspect([1 1 1]); axis off;
    view(3); view(viewangle); 
    hl=camlight;
    hl=camlight(hl,'headlight'); 
    lighting gouraud;
