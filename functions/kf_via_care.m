function [cntrl,P,L] = kf_via_care(sys,varargin)
% Lina notation
% no D-iggity 

    if isempty(varargin)
        
        A = sys.A;
        C = sys.C;
        W = sys.W;
        V = sys.V;
        
    else
        A = sys;
        C = varargin{1};
        W = varargin{2};
        V = varargin{3};
        
    end
    
    [ny,nx] = size(C);
    
% compute control and observer gains    
    
% observer
    P = icare(A',C',W,V,0,eye(nx),0);
    L = P*C'/V;    
    
% assemble dynamic controller    
    Ak = A - L*C;
    Bk = L;
    Ck = C;
    Dk = zeros(ny,ny);

    cntrl = ss(Ak,Bk,Ck,Dk);    
    

end