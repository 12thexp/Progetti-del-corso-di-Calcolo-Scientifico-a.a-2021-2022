function sol = thomas(A,b)
    % Solves linear algebraic equation where the coefficient matrix is 
    % tridiagonal. ld, md and ud stands for lower-, main, and upper-
    % diagonal respectively. a is the answer matrix and h is the solution.
    ld = diag(A,-1); md = diag(A); ud = diag(A,1);
    N = length(md) ;
    w = zeros(N, 1) ; g = zeros(N, 1) ;
    w(1) = ud(1)/md(1) ; g(1) = b(1)/md(1) ;
    if isrow(ud)
        ud = ud' ;
    end
    if isrow(ld)
        ld = ld' ;
    end
    ud = [ud; 0] ; ld = [0; ld] ;
    for i=2:N
        w(i) = ud(i)/(md(i)-ld(i)*w(i-1)) ;
        g(i) = (b(i)-ld(i)*g(i-1))/(md(i)-ld(i)*w(i-1)) ;
    end
    sol = zeros(N, 1) ;
    sol(N) = g(N) ;
    for i=N-1:-1:1
        sol(i) = -w(i)*sol(i+1)+g(i) ;
    end
end