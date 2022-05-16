function J = grad_extended_Phi_map_fd(Ac,Bc,Cc,Sigma,sys,varargin)


    [nc,nfi] = size(Bc);
    nfo = size(Cc,1);
    
    Sigma_inv = inv(Sigma);
    
    X = Sigma_inv(1:nc,1:nc);
    Y = Sigma(1:nc,1:nc);
    Ybar = Sigma(nc+1:end,nc+1:end);
    
    U = Sigma_inv(1:nc,nc+1:end);
    V = Sigma(1:nc,nc+1:end);
    
% nominal values
    [~,~,~,~,~,v0] = extended_Phi_map(Ac,Bc,Cc,Sigma,sys);

    J = [];
    
    h = 1e-6;
    
% w.r.t. A
    for i = 1:nc*nc
       
        Atmp = Ac;
        Atmp(i) = Atmp(i) + h;

        [~,~,~,~,~,v] = extended_Phi_map(Atmp,Bc,Cc,Sigma,sys);
                
        J = [J, (v - v0)/h];
                
    end
    
% w.r.t. B
    for i = 1:nc*nfi
       
        Btmp = Bc;
        Btmp(i) = Btmp(i) + h;

        [~,~,~,~,~,v] = extended_Phi_map(Ac,Btmp,Cc,Sigma,sys);
        
        J = [J, (v - v0)/h];
                
    end  
    
% w.r.t. C
    for i = 1:nfo*nc
       
        Ctmp = Cc;
        Ctmp(i) = Ctmp(i) + h;

        [~,~,~,~,~,v] = extended_Phi_map(Ac,Bc,Ctmp,Sigma,sys);
        
        J = [J, (v - v0)/h];
                
    end      
    
% w.r.t. Sigma - split it up

% w.r.t Y 
    
    for j = 1:nc

        for i = j:nc 

            Stmp = Y;

            if i == j

                Stmp(i,i) = Stmp(i,i) + h;

            else

                Stmp(i,j) = Stmp(i,j) + h;
                Stmp(j,i) = Stmp(j,i) + h;

            end           

            Simga_tmp = [Stmp, V; V', Ybar];
            
            [~,~,~,~,~,v] = extended_Phi_map(Ac,Bc,Cc,Simga_tmp,sys);

            J = [J, (v - v0)/h];

        end

    end 
       
  
% w.r.t V 
    
    for j = 1:nc

        for i = 1:nc 

            Stmp = V;

            Stmp(i,j) = Stmp(i,j) + h;
           
            Simga_tmp = [Y, Stmp; Stmp', Ybar];
            
            [~,~,~,~,~,v] = extended_Phi_map(Ac,Bc,Cc,Simga_tmp,sys);

            J = [J, (v - v0)/h];

        end

    end      
    
    
% w.r.t Ybar 
    
    for j = 1:nc

        for i = j:nc 

            Stmp = Ybar;

            if i == j

                Stmp(i,i) = Stmp(i,i) + h;

            else

                Stmp(i,j) = Stmp(i,j) + h;
                Stmp(j,i) = Stmp(j,i) + h;

            end           

            Simga_tmp = [Y, V; V', Stmp];
            
            [~,~,~,~,~,v] = extended_Phi_map(Ac,Bc,Cc,Simga_tmp,sys);

            J = [J, (v - v0)/h];

        end
        
    end

   
    
end

