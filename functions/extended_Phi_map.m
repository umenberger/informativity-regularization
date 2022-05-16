function [K,L,M,X,Y,v] = extended_Phi_map(Ac,Bc,Cc,Sigma,sys)

    A = sys.A;
    B = sys.B;
    C = sys.C;

    nx = size(A,1);
    ny = size(C,1);

    nc = size(Ac,1);
    nci = size(Bc,2);
    nco = size(Cc,1);

    % now switch to SCHERER NOTATION, don't mixup with Lina notation
    % - Sigma is inv(cal X)

    Y = Sigma(1:nx,1:nx);
    V = Sigma(1:nx,nx+1:end);

    calX = inv(Sigma);
    X = calX(1:nx,1:nx);
    U = calX(1:nx,nx+1:end);

    KLMN = [U X*B; zeros(nx,nco) eye(nco)]*[Ac Bc; Cc zeros(nco,nci)]*...
           [V' zeros(nc,ny); C*Y eye(ny)] + ...
           blkdiag(X*A*Y,zeros(nco,nci));

    K = KLMN(1:nc,1:nc);
    L = KLMN(1:nc,nc+1:end);
    M = KLMN(nc+1:end,1:nc);
%     N = KLMN(nc+1:end,nc+1:end);
    
    v = [K(:); L(:); M(:); sym2vec(X); sym2vec(Y)];

end


