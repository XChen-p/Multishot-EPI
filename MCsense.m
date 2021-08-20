function im = MCsense(k_obj, sample,sense,par)
para.sample=sample;
para.sens=sense;
[para.kx,para.ky,para.coil]=size(k_obj);

    % x subproblem
    % argmin_x { (0.5)||Mx-k||_2^2  }
    %     % (M'M x = M'k 
 
    
    para.y   =   reshape((1/2)*M_adj(k_obj,para),[],1);
    %pcg
    [im,~]   =   pcg(@(x)cgfun(x, para), para.y, 1E-6, 1000);
%      k      =   fftdim(reshape(im, para.kx, para.ky),1:2);
   im       =   abs(reshape(im, para.kx, para.ky));
     
    
end

function k_mc = M_fwd(im_sc, para)
% S*F*Sens
 
 
im_mc = repmat(im_sc,[1,1,para.coil,1]).*para.sens;
k_mc  = fftdim(im_mc,1:2).*para.sample;

% k_mc=(F*im_mc).*para.sample;

end
function im_sc = M_adj(k_mc, para)
%  sens*F*S
% k_mc=k_mc.*para.sample;
% for i=1:par.coil
%     im_mc(:,:,i)=par.FT'*k_mc(:,:,i);
% end

im_mc = ifftdim(k_mc.*para.sample,1:2);%n*l*coil*shot
im_sc = sum(im_mc.*conj(para.sens),3);%n*l*1*shot

end

function q = cgfun(x, para)
% S*F*Sens*
x = reshape(x, para.kx, para.ky); 
q = 1/2*M_adj(M_fwd(x,para),para) ;
q = reshape(q,[],1);

 
end