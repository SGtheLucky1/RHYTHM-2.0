% octree basis

% These are unit verticies
octb=[0,0,1/2
     0,0,0
     1/2,0,0
     1/2,0,1/2
     0,1/2,1/2
     0,1/2,0
     1/2,1/2,0
     1/2,1/2,1/2
     1,0,1/2
     1,0,0
     1,1/2,0
     1,1/2,1/2
     0,0,1
     1/2,0,1
     1/2,1/2,1
     0,1/2,1
     1,0,1
     1,1/2,1
     0,1,0
     1/2,1,0
     0,1,1/2
     1/2,1,1/2
     1,1,1/2
     1,1,0
     0,1,1
     1/2,1,1
     1,1,1];
     
% These are which verticies belong to each voxel, numbered
% counter-clockwise starting with from, lower, left vertex
octbi=[3,0,1,0,0,0,0,0
       1,0,0,0,0,0,0,0
       2,1,0,0,0,0,0,0
       4,3,2,1,0,0,0,0
       7,0,5,0,3,0,1,0
       5,0,0,0,1,0,0,0
       6,5,0,0,2,1,0,0
       8,7,6,5,4,3,2,1
       0,4,0,2,0,0,0,0
       0,2,0,0,0,0,0,0
       0,6,0,0,0,2,0,0
       0,8,0,6,0,4,0,2
       0,0,3,0,0,0,0,0
       0,0,4,3,0,0,0,0
       0,0,8,7,0,0,4,3
       0,0,7,0,0,0,3,0
       0,0,0,4,0,0,0,0
       0,0,0,8,0,0,0,4
       0,0,0,0,5,0,0,0
       0,0,0,0,6,5,0,0
       0,0,0,0,7,0,5,0
       0,0,0,0,8,7,6,5
       0,0,0,0,0,8,0,6
       0,0,0,0,0,6,0,0
       0,0,0,0,0,0,7,0
       0,0,0,0,0,0,8,7
       0,0,0,0,0,0,0,8];
      