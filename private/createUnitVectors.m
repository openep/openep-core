function uv=createUnitVectors(v)

m=sqrt(v(:,1).^2+v(:,2).^2+v(:,3).^2);
uv=[v(:,1)./m,v(:,2)./m,v(:,3)./m];

end