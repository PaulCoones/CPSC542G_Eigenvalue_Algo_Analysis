function [lambda, itn] = qreig(A, tol, max_iter)
  %
  % function [lambda,itn] = qreig(A,Tol)
  %
  % Find all real eigenvalues lambda of A
  % Return also iteration counters in itn

  arguments
      A
      tol = 1e-8
      max_iter = 1000
  end

  % First stage, bring to upper Hessenberg form
  % A = houseeig(A); % NOT NEEDED FOR TriDiag Matrices
  

  % Second stage: deflation loop
  n = size(A,1);
  lambda = [];
  itn = [];

  for j = n:-1:1
    % Find jth eigenvalue
    [lambda(j),itn(j),A] = qrshift(A(1:j,1:j),tol);
  end

  function [lam,iter,A] = qrshift(A,tol)
    %
    % function [lam,iter,A] = qrshift(A,tol)
    %
    % Find one eigenvalue lam of A in upper Hessenberg form,
    % return iteration count, too. Also improve A for future

    m = size(A,1);
    lam = A(m,m);
    iter = 0;
    I = eye(m);

    if m == 1
        return
    end

    while (iter < max_iter) % Max number of iterations

        if (abs(A(m,m-1)) < tol)
            return
        end

        iter = iter+1;

        [Q,R] = qr(A-lam*I); % Compute the QR decomposition
        A = R*Q+lam*I;       % Find the next iterate

        lam = A(m,m);       % Next shift
    end
  end
end
