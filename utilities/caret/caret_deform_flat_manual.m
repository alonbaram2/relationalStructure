function [C]=caret_deform_flat_manual(C,T,varargin);
% Aligns the sphere of the 'white matter surface' of the suit atlas
% VARARGIN:
%   Each varargin is a step
%   {'manual'}: Manual Gaussian shift, until empty
%   {'gaussian',start,end,radius}: Gaussian shift
%   {'latitude',start,end,radius}: distort latitude
%   {'vertex',v1,v2}: Attempts to move v1 on v2

% Make a topo structure with each Node as a first one in the triangle
numvert=size(C,1);
fprintf('Gettinging tiles......');
for i=1:numvert
    Tiles{i}=[];
    for j=1:3
        Tiles{i}=[Tiles{i};find(T(:,j)==i)];
    end;
end;


done=0;
n=1;

while ~done
    
    % Plot current state
    hold off;
    plot(C(:,1),C(:,2),'k.'); hold on;
    for i=1:3
        i1=i;
        i2=mod(i,3)+1;
        LX=[C(T(:,i1),1) C(T(:,i2),1)];
        LY=[C(T(:,i1),2) C(T(:,i2),2)];
        line(LX',LY','Color',[0 0 0]);
    end;
    set(gca,'XLim',[-2.5 -2],'YLim',[-1 -0.4]);
    
    % Area
    V1=C(T(:,3),:)-C(T(:,1),:);
    V2=C(T(:,2),:)-C(T(:,1),:);
    V3=cross(V1,V2);
    A=sqrt(sum(V3.*V3,2))/2;  % Area of the triangle
    A=V3(:,3);
    
    
    if n>length(varargin)
        done=1;
        break;
    end;
    
    % Derive command for the next step
    command=varargin{n};
    switch(command{1})
        case 'manual'
            % Get start and end
             x=ginput(2);   % Get direction of arrow
             if (isempty(x));
                 done=1;
                 break;
            end;
            % start = [-2.2601;-0.8207;0];
            % goal =  [-2.2547;-0.7270;0];
            start=[x(1,:) 0]';
            goal=[x(2,:) 0]';
            vect=goal-start;
            quiver(start(1),start(2),vect(1),vect(2),0);
            r=[ginput(1) 0]';  % Next click is the radius
            radius=norm(r-start);
            radius=0.12;
            d=surfing_eucldist(start,C')';
            w=exp(-1/(2*radius^2)*d.^2);
        case 'gaussian'
            start=[command{2} 0]';
            goal=[command{3} 0]';
            radius=command{4};
            vect=goal-start;
            quiver(start(1),start(2),vect(1),vect(2),0);
            d=surfing_eucldist(start,Sp1')';
            w=exp(-1/(2*radius^2)*d.^2);
            n=n+1;
        case 'latitude'
            start=command{2};
            goal=command{3};
            radius=command{4};
            vect=goal-start;;
            vect(1)=0; % Not in this direction
            d=Sp1(:,2)-start(2);
            w=exp(-1/(2*radius^2)*d.^2);
            n=n+1;
        case 'vertex'
            v1=command{2};
            v2=command{3};
            radius=command{4};
            start=[Sp1(v1,1:2) 0]';
            goal=[Sp2(v2,1:2) 0]';
            vect=goal-start;
            quiver(start(1),start(2),vect(1),vect(2),0);
            d=surfing_eucldist(start,Sp1')';
            w=exp(-1/(2*radius^2)*d.^2);
            n=n+1;
            keyboard;
    end;
    
    Cn=C+w*vect';
    C=caret_fixsurface_flat_old(Cn,T,Tiles,'fig',0);
    
    
    reduce=zeros(numvert,1);
    w=w.*(1-reduce*0.1);
    keyboard;
end;


