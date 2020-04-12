function badTriangles = caret_cecksphere( S1,T1)
% Check consistent alignement of triangles
%  function badTriangles = caret_cecksphere( S1,T1)

[S1n(:,1),S1n(:,2),S1n(:,3)]=sph2cart(Sp1(:,1),Sp1(:,2),Sp1(:,3));
V1=S1n(T1(:,2),:)-S1n(T1(:,1),:);
V2=S1n(T1(:,3),:)-S1n(T1(:,1),:);
V3=cross(V1,V2);
V3=bsxfun(@rdivide,V3,sqrt(sum(V3.*V3,2)));
S1nPl=S1n(T1(:,1),:)+V3;
r=sqrt(sum(S1nPl.*S1nPl,2));
find(r<100);
