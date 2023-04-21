function W = AssignGridWeights_8thNeighbour(img)
% A function that assigns grid weights in 8-Neighbour strategy, as discussed in our paper

[N1, N1, ~] = size(img);

N = N1^2;
mat = reshape(vec(1:N), [sqrt(N),sqrt(N)])'; %[1 2 3; 4 5 6; 7 8 9];                 % Sample matrix
[r, c] = size(mat);                          % Get the matrix size
diagVec1 = repmat([ones(c-1, 1); 0], r, 1);  % Make the first diagonal vector
                                             %   (for horizontal connections)
diagVec1 = diagVec1(1:end-1);                % Remove the last value
diagVec2 = [0; diagVec1(1:(c*(r-1)))];       % Make the second diagonal vector
                                             %   (for anti-diagonal connections)
diagVec3 = ones(c*(r-1), 1);                 % Make the third diagonal vector
                                             %   (for vertical connections)
diagVec4 = diagVec2(2:end-1);                % Make the fourth diagonal vector
                                             %   (for diagonal connections)
adj = diag(diagVec1, 1)+...                  % Add the diagonals to a zero matrix
      diag(diagVec2, c-1)+...
      diag(diagVec3, c)+...
      diag(diagVec4, c+1);
adj = adj+adj.';                             % Add the matrix to a transposed copy of
                                             %   itself to make it symmetric
G = gsp_2dgrid(sqrt(N));
Coord = G.coords; 
G = gsp_graph(adj);
G.coords = [(Coord(:,2)), flip(Coord(:,1))]; 
W_grid = G.W;

W = W_grid;

[N, N, ~] = size(img);

for n1 = 1 : N^2 - 1
    
    j1 = ceil(n1/N);
    
    i1 = n1 - (j1 - 1) * N;
    
    for n2 = n1 + 1 : N^2
        
        j2 = ceil(n2/N);

        i2 = n2 - (j2 - 1) * N;
       
        if W_grid(n1, n2) ~= 0
                       
            W(n1, n2) = exp(-norm(vec(img(i1, j1, :))-vec(img(i2, j2, :)), 'fro')^2);
        end
        
    end
    
end

W = W + W';

W_grid = [];

end