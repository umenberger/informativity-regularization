function [f,df,Z] = tr_inv_ct_kz(cntrl,sys)

% we will use Lina's ordering

[nx,nu] = size(sys.B);
ny = size(sys.C,1);

nc = nx; % should make this more general...

% input the controller as a matrix 
Ac = cntrl(1:nc,1:nc);
Bc = cntrl(1:nc,nc+1:end);
Cc = cntrl(nc+1:end,1:nc);
Dc = cntrl(nc+1:end,nc+1:end);

A = sys.A;
B = sys.B;
C = sys.C;
W = sys.W;
V = sys.V;

% nb: B should be zero for filtering applications
Acl = [A B*Cc; Bc*C Ac]; %Acl

% df = zeros(nu+nc,ny+nc);

% Check for closed-loop stability      
    if max(real(eig(Acl))) < 0
        
        % Closed-loop ctrb Gramian
        MWV = [W zeros(nx,nc); zeros(nc,nx) Bc*V*Bc'];  
        X = lyap(Acl,MWV); % Sigma
        X11 = X(1:nx,1:nx); % Sigma_x
        X12 = X(1:nx,nx+1:end); 
        X21 = X12'; % Sigma_xhatx
        X22 = X(nx+1:end,nx+1:end);

        Z = (X12/X22)*X12';
        
% cost function        
        f = trace(inv(Z));
        
% gradient
        X12inv = eye(nc)/X12;
        kew = lyap(Acl',blkdiag(zeros(nx),X12inv/X21));
        
        oh_mat = [zeros(nx), (X21\X22)*(X12inv/X21);
                  (X12inv/X21)*(X22/X12), zeros(nc)];
        
        oh = lyap(Acl',oh_mat);
        
        df = 2*[zeros(nx), eye(nx)]*(kew - oh)*[X12, X11*C';
                                                X22, X21*C' + Bc*V];
        
    else
        
        f = inf;
        df = nan;
        
    end
                                            
end

