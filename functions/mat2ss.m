function sys = mat2ss(M,n)

    sys = ss(M(1:n,1:n),M(1:n,n+1:end),M(n+1:end,1:n),M(n+1:end,n+1:end));

end