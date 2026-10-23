function A = tridiag(n, L)
    h = L/(n+1);
    A = (2*eye(n) - diag(ones(n-1,1),1) - diag(ones(n-1,1),-1)) / h^2;
end