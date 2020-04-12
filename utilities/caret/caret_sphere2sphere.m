function P=caret_sphere2sphere(V1,F1,V2);
% Resamples one sphere into the other
% Returns barycentric coordinates in sphere 1 for each of the nodes insphere 2
% INPUT:
% V1: Sphere 1 coordinates of vertices
% F1: Faces of sphere 1
% V2: Sphere 2 coordinates of vertices
% OUTPUT
% P.indx: Nodes in sphere1
% P.lambda: Lambda weights on nodes of sphere1
n2f=surfing_nodeidxs2faceidxs(F1');
strtThres = 0.12;

% Preallocat answer structure for speed
numV=size(V2,1);
P.indx=zeros(numV,3);
P.lambda=zeros(numV,3);

% Loop over nodes of surface 2 and figure out corresponding
% coordinates on Surface 1
for v=1:numV
    thresh = strtThres;
    
    if (mod(v,100)==0)
        fprintf('.'); 
    end; 
    
    % Distance from every node
    lambdas=surfing_eucldist(V2(v,:)',V1');
    R = sqrt(sum(V2(v,:).^2)); % Radius of the sphere
    
    closestFace=[];
    closestIndex = 0;
    found = 0;
    
    while ~found
        if thresh > R*2 % Know when to give up (should never happen)
            theFace = NaN;
            lambda = NaN;
            warning('This is bad.');
            keyboard;
        end
        % Get a few of the closest vertices
        possibleVertices = lambdas<R*thresh;
        
        possibleFaces = unique(n2f(possibleVertices,:));
        possibleFaces = possibleFaces(possibleFaces>0);
        % avgFaceDist = sum(lambdas(F1(possibleFaces,:)),2);
        % [~, closestFace] = sort(avgFaceDist);
        thresh = thresh*5;
        for f=1:length(possibleFaces)
            theFace = possibleFaces(f);
            fV = V1(F1(theFace,:),:);
            
            % You can divide the sphere into 3 half-spaces defined by the origin of
            % the sphere and 2 of the 3 vertices of the triangle. If the direction
            % from the half-spaces to the triangle has the same sign, then the
            % point lies inside the triangle. Apparently, the determinant of a
            % matrix can give you this info.
            %     http://forum.beyond3d.com/archive/index.php/t-48658.html
            d1 = det([fV(1:2,:); V2(v,:)]);
            d2 = det([fV(1,:); V2(v,:); fV(3,:)]);
            d3 = det([V2(v,:); fV(2:3,:)]);
            ds = [d1 d2 d3];
            if (all(ds>0) || all(ds<0))
                found=1;
                % Calculate BArycentric coordinates:
                x=[V2(v,:)-fV(1,:)]';  % Vector to projects
                T=[fV(2,:)-fV(1,:);fV(3,:)-fV(1,:)]';
                lmbd=[pinv(T)*x]';
                break;
            end;
        end;
    end;
    P.indx(v,:) = F1(theFace,:);
    P.lambda(v,:) = [1-sum(lmbd) lmbd];
end;
fprintf('done\n'); 

