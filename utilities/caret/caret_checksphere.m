function [badTriangles,Tfixed] = caret_checksphere( S1,T1)
% Check consistent alignement of triangles
%  function badTriangles = caret_cecksphere( S1,T1)
radius=sqrt(S1(1,:)*S1(1,:)');  % Determine constant radius 
V1=S1(T1(:,2),:)-S1(T1(:,1),:);
V2=S1(T1(:,3),:)-S1(T1(:,1),:);
V3=cross(V1,V2);
V3=bsxfun(@rdivide,V3,sqrt(sum(V3.*V3,2)));
S1nPl=S1(T1(:,1),:)+V3;
r=sqrt(sum(S1nPl.*S1nPl,2));
badTriangles=find(r<radius); % Inward pointing vectors
Tfixed=T1; 
if (length(badTriangles)>0)
    d=Tfixed(badTriangles,2); 
    Tfixed(badTriangles,2)=Tfixed(badTriangles,3); 
    Tfixed(badTriangles,3)=d; 
end; 