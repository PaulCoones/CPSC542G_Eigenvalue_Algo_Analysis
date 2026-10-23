% Use to directly compare two algorithms on specific sets of n and L
% values

clear all;
close all;

% Inits
n1=25;
n2=50;
n3=75;
n4=100;
n5=500;
l1 = 1;
l2 = 1;
l3 = 1;
l4 = 1;
l5 = 1;

% Matrices
Z1 = triDiag(n1,l1);
Z2 = triDiag(n2,l2);
Z3 = triDiag(n3,l3);
Z4 = triDiag(n4,l4);
Z5 = triDiag(n5,l5);

% Init relevant eigenvalue vectors
h1 = spacingVec(n1);
h2 = spacingVec(n2);
h3 = spacingVec(n3);
h4 = spacingVec(n4);
h5 = spacingVec(n5);
k1 = eigenStorage(n1);
k2 = eigenStorage(n2);
k3 = eigenStorage(n3);
k4 = eigenStorage(n4);
k5 = eigenStorage(n5);

% Compute exact analytical eigenvalues
lambda1 = eigenTriExact(h1, k1, n1);
lambda2 = eigenTriExact(h2, k2, n2);
lambda3 = eigenTriExact(h3, k3, n3);
lambda4 = eigenTriExact(h4, k4, n4);
lambda5 = eigenTriExact(h5, k5, n5);

% Compare for each pair of (n, L) values
lambdaBi1 = eigenBisection(Z1); % Change Bisection algo to X algo
lambdaJ1 = eigenJacobi(Z1);     % Change Jacobi algo to Y algo
errB1 = norm(lambda1 - lambdaBi1, 2)
errJ1 = norm(lambda1 - lambdaJ1, 2)

lambdaBi2 = eigenBisection(Z2);
lambdaJ2 = eigenJacobi(Z2);
errB2 = norm(lambda2 - lambdaBi2, 2)
errJ2 = norm(lambda2 - lambdaJ2, 2)

lambdaBi3 = eigenBisection(Z3);
lambdaJ3 = eigenJacobi(Z3);
errB3 = norm(lambda3 - lambdaBi3, 2)
errJ3 = norm(lambda3 - lambdaJ3, 2)

lambdaBi4 = eigenBisection(Z4);
lambdaJ4 = eigenJacobi(Z4);
errB4 = norm(lambda4 - lambdaBi4, 2)
errJ4 = norm(lambda4 - lambdaJ4, 2)

lambdaBi5 = eigenBisection(Z5);
lambdaJ5 = eigenJacobi(Z5);
errB5 = norm(lambda5 - lambdaBi5, 2)
errJ5 = norm(lambda5 - lambdaJ5, 2)

function h = spacingVec(n)
    h = 1/(n+1);
end

function k = eigenStorage(n)
    k = (1:n)';
end

function lambda = eigenTriExact(h, k, n) 
    lambda = (4/h^2) * sin(k*pi/(2*(n+1))).^2;
end