A1 = [0.6 -0.2 0.3 -0.8; 1 0 0 0; 0 0 0 0; 0 0 1 0];
B1 = [0; 0; 1; 0];

% Diagonalizable?
[V, D] = eig(A1);
rV = rank(V);
if rV == size(A1,1)
    disp('Dynamic matrix is diagonalizable.')
else
    disp('Dynamic matrix is not diagonalizable.')
end