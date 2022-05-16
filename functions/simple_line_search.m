function [t, xnew, problem] = simple_line_search(x, f, g, d, gamma, c, func)

% x - current point
% f - function value at x
% g - gradient at x
% d - search direction 
% gamma, c - line search parameters
% f - function to be optimized

gtd = trace(g'*d);

t = 1;
% t = 1e-1;

xnew = x + t*d;

fnew = func(xnew);

problem = 0;

    while fnew - f >= t*c*gtd

        t = gamma*t;

        xnew = x + t*d;
        
        fnew = func(xnew);
        
        if t < 1e-15
            fprintf('\tLine search warning: t = %.3e\n',t)
            problem = 1;
            break
        end

    end
    

end