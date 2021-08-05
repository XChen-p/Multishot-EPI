function im = ADMM(k_obj, x0, par)
% 
% Initialise x, z, y
x = x0;

z = Hankel(x, par.f, '2D');
y = z*1;
x0= x;
Hx= z;
g = z;

% H'H is just voxel-wise scaling by N
% i.e., (H'*H)*x = N.*x
[~, par.N]  =   pinv_hankel(Hankel(ones(par.kx,par.ky,par.coilcom,par.shot), par.f, '2D'), par.f, par.kx, par.ky, par.coilcom,par.shot);

% Define shortcut functions for hankel and adjoint
H_fwd   =   @(x)Hankel(x, par.f, '2D');
H_adj   =   @(x)pinv_hankel(x, par.f, par.kx, par.ky, par.coilcom, par.shot).*par.N;

for i = 1:par.niter
    
    % z subproblem
    % argmin_z { lambda||z||_* + (rho/2)||z - H(x) + y||_2^2
    % Solve via singular value soft thresholding
    [u,s,v] =   svd(Hx - y, 'econ');
%     s(s<par.lambda/par.rho)=0; % hard thresolding truncation
%     s      =   diag(max(diag(s) - par.lambda/par.rho, 0)); % soft thresolding shrinkage

    z=u*s*v';
    
    z_hat   =   par.r*z + (1-par.r)*g;
    
        
%     x subproblem
%     argmin_x { (1/2)||Mx-k||_2^2 + (rho/2)||z - H(x) + y||_2^2}
        % (M'M + (rho)H'H)x = M'k + (rho)H'(z+y)

    par.y   =   reshape((1/2)*M_adj(k_obj,par)+ (par.rho/2)*H_adj(z_hat+y) ,[],1); 
% pcg
    [k,~]   =   pcg(@(x)cgfun(x, par), par.y, 1E-6, 100); 
    x       =   reshape(k, par.kx, par.ky, par.coilcom,[]);
    
    Hx      =   H_fwd(x);    

    
    % dual update
    y       =   y + z_hat - Hx;
 
    
    % penalty adjustment
    s  = par.rho*H_fwd(x-x0);
    if norm(z(:)-Hx(:)) > 10*norm(s(:))
        par.rho = par.rho*par.m;
        y       = y/par.m;
    elseif norm(s(:)) > 10*norm(z(:)-Hx(:))
        par.rho = par.rho/par.m;
        y       = y*par.m;
    end
    
    
    x0 = x;
    g  = Hx;
    
    
    % print cost
    if par.verbose
        cost(i, x, Hx, k_obj, par);      
    end
    
end

im  = sos(ifftdim(reshape(x,[par.kx,par.ky,par.coilcom*par.shot]),1:2))./sqrt(par.state);
end

function k_mc = M_fwd(k_sc, par)
im_sc = ifftdim(k_sc,1:2);%n*l*1*shot
im_mc = repmat(im_sc,[1,1,par.coil,1]).*par.sens;
k_mc  = fftdim(im_mc,1:2).*par.sample;
end

function k_sc = M_adj(k_mc, par)
im_mc = ifftdim(k_mc.*par.sample,1:2);%n*l*coil*shot
im_sc = sum(im_mc.*conj(par.sens),3);%n*l*1*shot
k_sc  = fftdim(im_sc,1:2);
end


function q = cgfun(x, par)

x = reshape(x, par.kx, par.ky, par.coilcom,[]);
q = (1/2)*M_adj(M_fwd(x,par),par) + (par.rho/2)*(par.N.*x); 
q = reshape(q,[],1);
end


function cost(i, x, Hx, k_obj, par)
[~,s,~] = svd(Hx, 0);

c1 = norm(reshape(M_fwd(x,par) - k_obj, [],1))^2;
c2 = sum(diag(s));

c = 1/2*c1+par.lambda*c2;
fprintf(1, 'Iter: %d  Cost: %f  c1: %f  c2: %f  Rho: %f\n', i, c, c1,c2,par.rho);
end
