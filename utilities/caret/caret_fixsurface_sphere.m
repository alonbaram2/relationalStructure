function C=caret_fixsurface_sphere(C,T,Tiles,varargin);
% Attempts to wiggle the triangles around a bad
% trian


fig=0;
vararginoptions(varargin,{'fig'});
radius=sqrt(sum(C(1,:).^2));



V1=C(T(:,3),:)-C(T(:,1),:);
V2=C(T(:,2),:)-C(T(:,1),:);
V3=cross(V1,V2,2);
Ctest=C(T1,:)+V3; 
r=sqrt(sum(Ctest.*Ctest,2)); 
badTriangles=(r<radius); 

if (isempty(badTriangles))
    return;
end;

% Set up the exploration steps
delta=[-1 0;1 0;0 -1;0 1];
X=[delta ones(4,1)];
P=inv(X'*X)*X';

% Only consider the triangles connection to the bad ones: speeds things up
TT=T(badTriangles,:);               % Bad triangles
iC=unique(TT(:));                   % Index of bad vertices
allTriangles=sort(unique(vertcat(Tiles{iC})));
TT=T(allTriangles,:);         % All triangles bordering bad triangles

% Now take steps in the right direction until fixed
done=0;
n=0;
numBad=length(badTriangles);
while ~done
    Cn=C;
    V1=C(TT(:,3),:)-C(TT(:,1),:);
    V2=C(TT(:,2),:)-C(TT(:,1),:);
    V3=cross(V1,V2,2);
    Ctest=C(T1,:)+V3; 
    r=sqrt(sum(Ctest.*Ctest,2)); 
    iB=find(r<radius); 

    iB=find(A<0);
    if (isempty(find(A<0)))
        break;
    else
        n=n+1;
        iB=(find(A<0));
        
        fprintf('Bad: %d %d\n',n,length(iB));
        if (numBad>length(iB))
            n=0;
            numBad=length(iB);
        elseif (n>50)
            hold off;
            plot(C(:,1),C(:,2),'k.'); hold on;
            % set(gca,'XLim',[-2.4 -2.1],'YLim',[-0.7 -0.6]);
            keyboard;
        end;
    end;
    deriv=zeros(size(Cn,1),2);
    [Sp(:,1),Sp(:,2),Sp(:,3)]=cart2sph(Cn(:,1),Cn(:,2),Cn(:,3));
    Spn=Sp; 
    for i=iC'
        t=Tiles{i};
        V1=Cn(T(t,3),:)-Cn(T(t,1),:);
        V2=Cn(T(t,2),:)-Cn(T(t,1),:);
        V3=cross(V1,V2,2);
        Ctest=C(T1,:)+V3; 
        r=sqrt(sum(Ctest.*Ctest,2)); 
        iB=find(r<radius); 
        Area=r-radius;
        if (min(Area)<0)
            explorestep=0.02; 
            for d=1:4
                Spn(i,1:2)=Spn(i,1:2)+delta(d,:)*explorestep;
                V1=Cn(T(t,3),:)-Cn(T(t,1),:);
                V2=Cn(T(t,2),:)-Cn(T(t,1),:);
                V3=cross(V1,V2,2);
                a=V3(:,3);
                a(a>0)=0;
                area(d,1)=mean(a);  % Area of the triangle
            end;
            d=P*area;
            deriv(i,1:2)=d(1:2);
            Cn(i,:)=C(i,:); % Reset the coordinate
        end;
    end;
    % deriv(iC,:)
    C(:,1:2)=C(:,1:2)+deriv;
end;