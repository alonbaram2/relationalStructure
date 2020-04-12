function P=freesurfer_sphere2sphere(sphere1,sphere2);
% Resamples one sphere into the other
% Returns barycentric coordinates in sphere 1 for each of the nodes in
% sphere 2
% OUTPUT
% P.indx: Nodes in sphere1
% P.lambda: Lambda weights on nodes of sphere1
fig=0;
[V1,F1]=read_surf(sphere1);
[V2,F2]=read_surf(sphere2);
F1=F1+1;
F2=F2+1;
n2f1=surfing_nodeidxs2faceidxs(F1');
n2f2=surfing_nodeidxs2faceidxs(F2');

% Loop over nodes of surface 2 and figure out corresponding
% coordinates on Surface 1
for i=1:size(V2,1)
    p=V2(i,:)';
    Dist=surfing_eucldist(p,V1');
    [DistS,I]=sort(Dist,2,'ascend'); 
    face=n2f1(I(1:6),:);
    face=unique(face(face>0)); 
    
    numFace=length(face);
    lmb=NaN(2,numFace);
    proj=NaN(3,numFace); 
    if (fig>0)
        plot3(p(1),p(2),p(3),'r.');hold on;
    end;
    for f=1:numFace
        FF=F1(face(f),:);
        x=[p-V1(FF(1),:)'];  % Vector to projects
        T=[V1(FF(2),:)-V1(FF(1),:);V1(FF(3),:)-V1(FF(1),:)]';
        lmb(:,f)=inv(T'*T)*T'*x;
        % proj(:,f)=T*lmb(:,f)+V1(FF(1),:)';
        % Visualization ]
        if (fig>0)
            C=V1(FF,:);
            for j=1:3
                plot3(C(j,1),C(j,2),C(j,3),'k.'); hold on;
            end;
            line(C(1:2,1)',C(1:2,2)',C(1:2,3)');
            line(C(2:3,1)',C(2:3,2)',C(2:3,3)');
            line(C([1 3],1)',C([1 3],2)',C([1 3],3)');
        end;
    end;
    goodF=find(sum(lmb)<1 & all(lmb>0));
    if (length(goodF)==0)
        warning('no triangle found .... taking closest one');
        II=find(all(lmb>0));
        if (isempty(II))
            goodF=1;
        else
            ll=sum(lmb(:,II));
            [~,j]=min(ll);
            goodF=II(j);
        end;
    end;
    if (length(goodF)>1)
        warning('two or more triangles found, taking first');
        goodF=goodF(1); 
    end;
    
    P.indx(i,:)=F1(face(goodF),:);
    P.lambda(i,:)=[1-sum(lmb(:,goodF)) lmb(:,goodF)'];

    % Draw in final projection 
    if (fig>0)
            X=V1(P.indx(i,:)',:);
            proj=sum(bsxfun(@times,X,P.lambda(i,:)'));

        
        plot3(proj(1),proj(2),proj(3),'b.');
        hold off; 
        axis equal;
        keyboard; 
    end;
    
    
    if (mod(i,100)==0)
        fprintf('%d\n',i);
    end;
end;

