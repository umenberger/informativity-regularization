function [f,df,X] = kf_state_cost_ct(cntrl,sys)

% In this function the cost is | x - Cc xhat |^2
% as opposed to | y - Cc xhat |^2, as in Kalman filtering.
% The input to the filter is still y, and C can be short/wide.

[ny,nx] = size(sys.C);

% This old calculation no longer works
% infer size of controller
% nc_ny = size(cntrl,2);
% nc = nc_ny - ny;
% Now the output dimension of the controller must be nx
nc_nx = size(cntrl,1);
nc = nc_nx - nx;

% input the controller as a matrix 
Ac = cntrl(1:nc,1:nc);
Bc = cntrl(1:nc,nc+1:end);
Cc = cntrl(nc+1:end,1:nc);
Dc = cntrl(nc+1:end,nc+1:end);

A = sys.A;
C = sys.C;
W = sys.W;
V = sys.V;

Acl = [A zeros(nx,nc); Bc*C Ac];

    if max(real(eig(Acl))) < 0

        BclBcl = [W zeros(nx,nc); zeros(nc,nx) Bc*V*Bc'];
        CclCcl = [eye(nx), -Cc; -Cc', Cc'*Cc];

        X = lyap(Acl,BclBcl);

        f = trace(CclCcl*X);
        
        if nargout > 1
            
            Y = lyap(Acl',CclCcl);   
            
        % decomp
            X11 = X(1:nx,1:nx);
            X12 = X(1:nx,nx+1:end);
            X22 = X(nx+1:end,nx+1:end);

%             Y11 = Y(1:nx,1:nx);
            Y12 = Y(1:nx,nx+1:end);
            Y22 = Y(nx+1:end,nx+1:end);

        % gradient
            dfdA = 2*(Y12'*X12 + Y22*X22);
            dfdB = 2*(Y22*Bc*V + Y22*X12'*C' + Y12'*X11*C');
            dfdC = 2*(-X12 + Cc*X22);
            dfdD = zeros(nx,ny);

            df = [dfdA, dfdB;
                  dfdC, dfdD];            
            
        end        
        
        
        
    else
        
        f = inf;
        df = nan;
        X = nan;
        
    end


end