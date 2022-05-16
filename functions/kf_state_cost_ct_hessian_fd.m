function H = kf_state_cost_ct_hessian_fd(cntrl,sys)

% dimensions
[ny,nx] = size(sys.C);
nc_nx = size(cntrl,1);
nc = nc_nx - nx;

% input the controller as a matrix 
Ac = cntrl(1:nc,1:nc);
Bc = cntrl(1:nc,nc+1:end);
Cc = cntrl(nc+1:end,1:nc);
Dc = cntrl(nc+1:end,nc+1:end);

nfi = size(Bc,2); % number of inputs to filter
nfo = size(Cc,1);

% w.r.t. Ac
H = [];

[~,df0,~] = kf_state_cost_ct(cntrl,sys);

% df0 = df0(1:nc,:); % we want w.r.t. Ck as well!

df0 = vectorize(df0);

h = 1e-5;

    for i = 1:nc*nc
       
        Atmp = Ac;
        Atmp(i) = Atmp(i) + h;
        
        cntrl_ = [Atmp, Bc; Cc, Dc];
        
        [~,df,~] = kf_state_cost_ct(cntrl_,sys);
        
        df = vectorize(df);
        
        H = [H, (df(:) - df0)/h];
           
    end
    
    for i = 1:nc*nfi
       
        Btmp = Bc;
        Btmp(i) = Btmp(i) + h;
        
        cntrl_ = [Ac, Btmp; Cc, Dc];
        
        [~,df,~] = kf_state_cost_ct(cntrl_,sys);
        
        df = vectorize(df);
        
        H = [H, (df(:) - df0)/h];
           
    end   
 
   
    for i = 1:nfo*nc
       
        Ctmp = Cc;
        Ctmp(i) = Ctmp(i) + h;
        
        cntrl_ = [Ac, Bc; Ctmp, Dc];
        
        [~,df,~] = kf_state_cost_ct(cntrl_,sys);
        
        df = vectorize(df);
        
        H = [H, (df(:) - df0)/h];
           
    end 
    
    
%     for i = 1:nfo*nfi
%        
%         Dtmp = Dc;
%         Dtmp(i) = Dtmp(i) + h;
%         
%         cntrl_ = [Ac, Bc; Cc, Dtmp];
%         
%         [~,df,~] = kf_state_cost_ct(cntrl_,sys);
%         
%         df = vectorize(df);
%         
%         H = [H, (df(:) - df0)/h];
%            
%     end     
    
    
    function v = vectorize(M)
       
% select components of gradient 
        MA = M(1:nc,1:nc);
        MB = M(1:nc,nc+1:end);
        MC = M(nc+1:end,1:nc);
        MD = M(nc+1:end,nc+1:end);
        
        v = [MA(:); MB(:); MC(:)];
%         v = [ MB(:); MA(:); MC(:)];

        
        
    end
    
    
    

end