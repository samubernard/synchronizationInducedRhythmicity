function [G,Left,Right]= numgrid3(R,n)
%NUMGRID Number the grid points in a three dimensional region.
%   G = NUMGRID3('R',n) numbers the points on an n-by-n grid in
%   the subregion of -1<=x<=1 and -1<=y<=1 and -1<=z<=1 determined by 'R'.
%   The regions available are:
%       'SCN'= 3D SCN, roughly two intersecting ellipsoids
%       'SCNslice'= 3D slice of the SCN along the VL-DM axis, 3 layers thick
%       '2Dslice'= 2D slice (1 layer thick)
%       'Neuron'= 1 single neuron (for compatibility)

G=zeros(n,n,n);
Left=zeros(n,n,n);
Right=zeros(n,n,n);
x = [-1, (-(n-3):2:(n-3))/(n-1), 1];
y = x;
z = y;
b = 0.45; % 0.45

switch R,
    case 'SCN'
        for i=1:n,
            for j=1:n,
                for k=1:n,
                    rightneuron = ((x(i)-b)^2/0.25 + y(j)^2 + z(k)^2/(.7)^2 < 1); % Right 
                    leftneuron = ((x(i)+b)^2/0.25 + y(j)^2 + z(k)^2 /(.7)^2< 1); % Left
                    G(i,j,k) = leftneuron || rightneuron;
                    Left(i,j,k) = leftneuron;
                    Right(i,j,k) = rightneuron;
                    
                end
            end
        end
    case 'SCNslice'
        for i=1:n,
            for j=1:n,
                for k=1:n,
                    rightneuron = ((x(i)-b)^2/0.25 + y(j)^2 + z(k)^2/(.7)^2 < 1); % Right 
                    leftneuron = ((x(i)+b)^2/0.25 + y(j)^2 + z(k)^2 /(.7)^2< 1); % Left
                    G(i,j,k) = leftneuron || rightneuron;
                    Left(i,j,k) = leftneuron;
                    Right(i,j,k) = rightneuron;
                    
                end
            end
        end
        S2=size(G,2);
        G=G(:,floor(S2/2)-1:floor(S2/2)+1,:);
        Left = Left(:,floor(S2/2)-1:floor(S2/2)+1,:);
        Right = Right(:,floor(S2/2)-1:floor(S2/2)+1,:);
    case '2Dslice'
        for i=1:n,
            for j=1:n,
                for k=1:n,
                    rightneuron = ((x(i)-b)^2/0.25 + y(j)^2 + z(k)^2/(.7)^2 < 1); % Right 
                    leftneuron = ((x(i)+b)^2/0.25 + y(j)^2 + z(k)^2 /(.7)^2< 1); % Left
                    G(i,j,k) = leftneuron || rightneuron;
                    Left(i,j,k) = leftneuron;
                    Right(i,j,k) = rightneuron;
                    
                end
            end
        end
        S2=size(G,2);
        G=G(:,floor(S2/2),:);
        Left = Left(:,floor(S2/2),:);
        Right = Right(:,floor(S2/2),:);
    case 'Neuron'
        G(1,1,1)=1;
        Left=G;
        Right=G;
end


k = find(G);
G = zeros(size(G));      % Convert from logical to double 3D array
G(k) = (1:length(k))';
k = find(Left);
Left = zeros(size(Left));      % Convert from logical to double 3D array
Left(k) = (1:length(k))';
k = find(Right);
Right = zeros(size(Right));      % Convert from logical to double 3D array
Right(k) = (1:length(k))';


