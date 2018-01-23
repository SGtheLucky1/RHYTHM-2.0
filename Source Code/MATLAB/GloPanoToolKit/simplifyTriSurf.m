%% Simplification Algorithm %%
clear all
cd('/Users/Chris/Data/PanoramicImaging/ExperimentsGWU/2017_0113_Rabbit/Optical')
load('06_projectedData')
faces = a.cells;
pts = a.pts;
centroids = a.centroids;
norms = a.norms;

% Create a variable of all the edges
edges = zeros(size(faces,1)*3,2);
for n = 1:size(faces,1)
   edges(3*(n-1)+1:3*(n-1)+3,:) = [faces(n,1) faces(n,2); 
       faces(n,1) faces(n,3); 
       faces(n,2) faces(n,3)];
end

% Sort the edges to match edges with different ordered vertices
edges = sort(edges,2);
% Establish counter
n = 1;
% Step through each edge finding and removing duplicates
while n < size(edges,1)
   % Create variable for checking for duplicates
   dupes = repmat(edges(n,:),[size(edges,1) 1]);
   % Check against current edges variable
   dupes = sum(dupes == edges,2) == 2;
   % If there are duplicates, remove all but the first one
   if sum(dupes)~=0
       % Grab the indices
       dupes = find(dupes);
       % Discard all but the first one
       dupes = dupes(2:end);
       edges(dupes,:) = [];
   end
   % Increment the counter
   n = n+1;
end

% Collect all normals and vertices adjacent to each edge
nbrCells = cell(size(edges,1),1);
nbrNorms = cell(size(edges,1),1);
nbrVerts = cell(size(edges,1),1);
ind = 1:size(edges,1);
[nbrCells,nbrNorms,nbrVerts] = findNeighborsStuff(ind,edges,faces,...
norms,nbrCells,nbrNorms,nbrVerts);
% % % for n = 1:size(edges,1)
% % %     % Find the neighboring cells
% % %     tmp = reshape(edges(n,:),[1 1 2]);
% % %     tmp = repmat(tmp,[size(faces,1) 3]);
% % %     nbrCells{n} = find(sum(sum(tmp == repmat(faces,[1 1 2]),2),3)~=0);
% % %     nbrNorms{n} = norms(nbrCells{n},:);
% % %     % Find the connected vertices
% % %     nbrVerts{n} = unique(faces(nbrCells{n},:));
% % %     tmp = repmat(nbrVerts{n},[1 2]) == repmat(edges(n,:),[size(nbrVerts{n},1) 1]);
% % % %     tmp = find(sum(tmp,2)~=0);
% % %     nbrVerts{n}(sum(tmp,2)~=0) = [];  
% % % end

% Calculate the minimizer (p), QEM[p], and QEM coefficients (a,b,c)
p = zeros(size(edges,1),3);
QEM = zeros(size(edges,1),1);
a = cell(size(edges,1),1);
b = cell(size(edges,1),1);
c = zeros(size(edges,1),1);
% Iterate through each edge and calculate the coefficients, p, and QEM
for n = 1:size(edges,1)
    % Coefficient a
    a{n} = [sum(nbrNorms{n}(:,1).^2) sum(nbrNorms{n}(:,1).*nbrNorms{n}(:,2)) sum(nbrNorms{n}(:,1).*nbrNorms{n}(:,3));
        sum(nbrNorms{n}(:,2).*nbrNorms{n}(:,1)) sum(nbrNorms{n}(:,2).^2) sum(nbrNorms{n}(:,2).*nbrNorms{n}(:,3));
        sum(nbrNorms{n}(:,3).*nbrNorms{n}(:,1)) sum(nbrNorms{n}(:,3).*nbrNorms{n}(:,2)) sum(nbrNorms{n}(:,3).^2)];
    % Coefficient b
    b{n} = [sum(nbrNorms{n}(:,1).*sum(nbrNorms{n}.*centroids(nbrCells{n},:),2));
        sum(nbrNorms{n}(:,2).*sum(nbrNorms{n}.*centroids(nbrCells{n},:),2));
        sum(nbrNorms{n}(:,3).*sum(nbrNorms{n}.*centroids(nbrCells{n},:),2))];
    % Coefficient c
    c(n) = sum(sum(nbrNorms{n}.*centroids(nbrCells{n},:),2).^2);
    % Point that minimizes QEM
    p(n,:) = transpose(a{n}\b{n});
    % Quadratic Error Metric (QEM)
    QEM(n) = p(n,:)*a{n}*transpose(p(n,:))-2*p(n,:)*b{n}+c(n);
end

% % % figure
% % % trisurf(faces,pts(:,1),pts(:,2),pts(:,3))
% % % h = waitbar(0,'Simplifying mesh...');
for m = 1:1000
    
    % 1) Identify the index of the incorrect neighborhood
    % 2) Track it to identify during which iteration it becomes incorrect
    
    
    disp(pts(5023,:))
    % Find the smallest QEM value
    [~,I] = min(QEM);
    rmEdge = edges(I,:);
    
% % %     % Visualize
% % %     figure
% % %     trisurf(faces,pts(:,1),pts(:,2),pts(:,3))
% % %     hold on
% % %     scatter3(pts(rmEdge,1),pts(rmEdge,2),pts(rmEdge,3),'ro','filled','SizeData',128)
% % %     scatter3(pts(nbrVerts{I},1),pts(nbrVerts{I},2),pts(nbrVerts{I},3),'go','filled','SizeData',128)
% % %     check1 = pts(nbrVerts{I},1);
% % %     check2 = pts(nbrVerts{I},2);
% % %     check3 = pts(nbrVerts{I},3);
% % %     hold off
    
    % Remove the edge
    edges(I,:) = [];
    coeffs{1} = a{I}; a(I) = [];
    coeffs{2} = b{I}; b(I) = [];
    coeffs{3} = c(I); c(I) = [];
    % Hold onto minimizer location, but remove minimizer and QEM from variables
    minimizer = p(I,:);
    pts = [pts; minimizer];
    p(I,:) = [];
    QEM(I) = [];
    nbrCells(I) = [];
    nbrNorms(I) = [];
    nbrVerts(I) = [];
    % Remove the triangles that share the edge (degenerate triangles) and their
    % associated norms
    rm = repmat(faces,[1 1 2]) == repmat(reshape(rmEdge,[1 1 2]),[size(faces,1) 3 1]);
    rm = find(sum(sum(rm,2),3) == 2);
    % Save out the triangle vertices to use for finding duplicate edges
    dupEdges = faces(rm,:);
    % Remove triangles and their norms
    faces(rm,:) = [];
    norms(rm,:) = [];
    centroids(rm,:) = [];
    % Update neighboring triangles with minimizer index
    tmp = reshape(rmEdge,[1 1 2]);
    tmp = repmat(tmp,[size(faces,1) 3]);
    nbrFaces = find(sum(sum(tmp == repmat(faces,[1 1 2]),2),3)~=0);
    col = (faces(nbrFaces,:) == rmEdge(1))+(faces(nbrFaces,:) == rmEdge(2));
    col = sum(col.*repmat(1:3,[size(col,1) 1]),2);
    ind = sub2ind([size(faces,1),size(faces,2)],nbrFaces,col);
    faces(ind) = size(pts,1);
    % Calculate new centroids at these faces
    centroids(nbrFaces,1) = sum(reshape(pts(faces(nbrFaces,:),1),[length(nbrFaces) 3]),2)/3;
    centroids(nbrFaces,2) = sum(reshape(pts(faces(nbrFaces,:),2),[length(nbrFaces) 3]),2)/3;
    centroids(nbrFaces,3) = sum(reshape(pts(faces(nbrFaces,:),3),[length(nbrFaces) 3]),2)/3;
    % Calculate new normals at these centroids
    v1 = pts(faces(nbrFaces,1),:)-centroids(nbrFaces,:);
    v2 = pts(faces(nbrFaces,2),:)-centroids(nbrFaces,:);
    for n = 1:size(v1,1)
        norms(nbrFaces(n),:) = cross(v2(n,:),v1(n,:));
        if sum(-centroids(nbrFaces(n),:).*norms(nbrFaces(n),:),2) > 0
            norms(nbrFaces(n),:) = cross(v1(n,:),v2(n,:));
        end
    end
    
    % Update edges with minimizer values
    nbrEdges = repmat(reshape(rmEdge,[1 1 2]),[size(edges,1) size(edges,2)]);
    nbrEdges = nbrEdges == repmat(edges,[1 1 2]);
    nbrEdges = find(sum(sum(nbrEdges,2),3));
    col = (edges(nbrEdges,:) == rmEdge(1))+(edges(nbrEdges,:) == rmEdge(2));
    col = sum(col.*repmat(1:2,[size(col,1) 1]),2);
    ind = sub2ind([size(edges,1),size(edges,2)],nbrEdges,col);
    edges(ind) = size(pts,1);
    
    % Remove redundant edges and associated QEM variables
    edges = sort(edges,2);
    % Establish a counter
    n = 1;
    % Create an empty variable for the duplicate indices
    dupes = [];
    while n < size(nbrEdges,1)
        % Create a variable for checking against the edges variable
        check = repmat(edges(nbrEdges(n),:),[size(nbrEdges,1) 1]);
        % Check against edges variable
        check = check == edges(nbrEdges,:);
        % Switch the self match to not a match (0)
        check(n,:) = 0;
        check = find(sum(check,2)==2);
        if ~isempty(check)
            % Grab the duplicates location
            dupes = [dupes;nbrEdges(check)];
            % Remove from the neighbor variable
            nbrEdges(check) = [];
        end
        % Increment counter
        n = n+1;
    end
    edges(dupes,:) = [];
    a(dupes) = [];
    b(dupes) = [];
    c(dupes) = [];
    p(dupes,:) = [];
    QEM(dupes) = [];
    nbrCells(dupes) = [];
    nbrNorms(dupes) = [];
    nbrVerts(dupes) = [];
    
    % Update QEM variables of remaining edges
    for n = 1:size(nbrEdges,1)
        % Grab relevant index
        ind = nbrEdges(n);
        % Coefficient a
        a{ind} = a{ind}+coeffs{1};
        % Coefficient b
        b{ind} = b{ind}+coeffs{2};
        % Coefficient c
        c(ind) = c(ind)+coeffs{3};
        % Point that minimizes QEM
        p(ind,:) = transpose(a{ind}\b{ind});
        % Quadratic Error Metric (QEM)
        QEM(ind) = p(ind,:)*a{ind}*transpose(p(ind,:))-2*p(ind,:)*b{ind}+c(ind);
    end
    
    % Update neighbors stuff
    [nbrCells,nbrNorms,nbrVerts] = findNeighborsStuff(nbrEdges,edges,faces,...
        norms,nbrCells,nbrNorms,nbrVerts);
% % %     
% % %     % Visualize
% % %     figure
% % %     trisurf(faces,pts(:,1),pts(:,2),pts(:,3))
% % %     hold on
% % %     scatter3(pts(end,1),pts(end,2),pts(end,3),'ro','filled','SizeData',128)
% % %     scatter3(check1,check2,check3,'go','filled','SizeData',128)
% % %     hold off
% % %     waitbar(m/1000)
end
close(h)

% Visualize
figure
trisurf(faces,pts(:,1),pts(:,2),pts(:,3))
