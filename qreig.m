function [lambda, itn] = qreig(A, tol, max_iter, rel_tol)
        % function [lambda,itn] = qreig(A,Tol)
        % Find all real eigenvalues lambda of A
        % Return also iteration counters in itn

        arguments
            A
            tol = 1e-8
            max_iter = 1000
            rel_tol = 0
        end
        
        % % Unneeded as A is tridiagonalized in our timings, and this
        % % pre tridiagonalization is timed aswell.
        % % First stage, bring to upper Hessenberg form
        % % IF it is not already in this form
        % if ~(isbanded(A, 1, size(A, 2)))
        %    fprintf("\n For Shifted QR Iteration - Matrix is not yet in upper Hessenburg form, reducing...\n")
        %    A = houseeig(A);
        % end

        % Second stage: deflation loop
        n = size(A,1);
        lambda = [];
        itn = [];
        for j = n:-1:1
            % Find jth eigenvalue
            [lambda(j),itn(j),A] = qrshift(A(1:j,1:j),tol,rel_tol);
        end

        function [lam,iter,A] = qrshift(A,tol,rel_tol)
            %
            % function [lam,iter,A] = qrshift(A,tol,rel_tol)
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

                % Determine whether to use a relative or absolute tolerance
                if rel_tol == 1
                    % Set the tolerance rel to neighbouring diag elem.
                    currentTol = tol * (abs(A(m-1,m-1)) + abs(A(m,m)));
                else
                    currentTol = tol;
                end

                % Check if the last eigenvalue has converged
                if (abs(A(m,m-1)) <= currentTol)
                    return
                end

                iter = iter+1;

                [Q,R] = qr(A-lam*I); % Compute the QR decomposition
                A = R*Q+lam*I;       % Find the next iterate

                lam = A(m,m);       % Next shift
            end
        end
    end