function cntrl_ = balance_X22(cntrl,sys)

nx = size(sys.A,1);

% [~,~,X] = kf_state_cost_ct(ss2mat(cntrl),sys);

if sum(sum(sys.B)) == 0 % state cost 
   
    [~,~,X] = kf_state_cost_ct(ss2mat(cntrl),sys);
    
else
    
    [~,~,X] = lqg_ct_lina(ss2mat(cntrl),sys);
    
end

try

L = chol(X(nx+1:end,nx+1:end));
L = L';

cntrl_ = cntrl;
cntrl_.A = L\cntrl.A*L;
cntrl_.B = L\cntrl.B;
cntrl_.C = cntrl.C*L;

catch
        
    Acl = [sys.A sys.B*cntrl.C; cntrl.B*sys.C cntrl.A];
    
    max(real(eig(Acl)))
    
    eig(Acl)
    
    debug = 1;
    
    
end

end