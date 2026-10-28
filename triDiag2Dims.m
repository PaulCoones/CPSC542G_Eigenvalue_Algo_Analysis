% Construct 2D finite-difference Laplacian
nx = 2^ceil(log2(n)/2);
ny = n/nx;

hx = L/(nx+1);
hy = L/(ny+1);

ex = ones(nx, 1);
ey = ones(ny, 1);

Kx = spdiags([-ex 2*ex -ex], -1:1, nx, nx) / hx^2;
Ky = spdiags([-ey 2*ey -ey], -1:1, ny, ny) / hy^2;

% Construct 2D Laplacian
A = kron(speye(ny), Kx) + kron(Ky, speye(nx));

% Tridiagonalize the matrix
Z = hess(full(A));