function [k_obj] = gen_obj(t,shot,noise_level)
load('ww_1_5.mat')
load('obj_1_5.mat')

ww=ww(:,(t-1)*shot+1:t*shot);
 
 
 [par.kx,par.ky,par.coil]=size(obj);

 
np=5;
p=ones(par.kx,par.ky,np);
for i=1:par.kx
        p(i,:,2)=(-49+i)/96;
end
for j=1:par.ky
        p(:,j,3)=(-49+j)/96;
end
for i=1:par.kx
    for j=1:par.ky
        p(i,j,4)=((-49+i)/96)^2-((-49+j)/96)^2;
    end
end
for i=1:par.kx
    for j=1:par.ky
        p(i,j,5)=(-49+i)/96*(-49+j)/96;
    end
end

phs=zeros(par.kx,par.ky,shot);
for i=1:shot
for j=1:np
    phs(:,:,i)=phs(:,:,i)+p(:,:,j)*ww(j,i);
end
end
phs=permute(repmat(phs,[1,1,1,par.coil]),[1,2,4,3]);


% Generate complex obj
obj = bsxfun(@times, obj, exp(1j*phs*1));
 

% generate k-space
k_obj   = fftdim(obj,1:2);

rng(t+100);
noise_i = randn(size(k_obj))*noise_level;
noise_j = randn(size(k_obj))*noise_level;
k_obj=noise_i+1j*noise_j+k_obj;

 

end
