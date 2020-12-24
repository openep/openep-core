function projU=projectVertexVectors(verts,u,normVerts)

%projects the velocities (u) defined on the vertices (verts) in the
%perpendicular to the vertex normal direction (normVerts):

projU=zeros(size(verts));
for i=1:length(verts)
    n=normVerts(i,:);
    ui=n*dot(n,u(i,:));
    projU(i,:)=u(i,:)-ui;
end

end