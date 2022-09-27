function [k_obj,par] = gen_obj(frame)
load('ww_1_5.mat')
load('obj_1_5.mat')
load('sens_1_5.mat')
% load('sens_s.mat')
% sens=padarray(sens,[2 0 0],'post');
% load('obj_phantom.mat') 
% load('sens_phantom.mat')
% sens=padarray(sens,[16 16 0],'both');
par.sens=sens; 
load('mask_11.mat')
SE=strel('sphere',3);
% mask=imdilate(mask,SE);
par.sens=sens;%.*repmat(mask,[1 1 32]);
ob=zeros(64,64);
ob(25:40,25:40)=1;
% par.mask=ob;
% xx=[1 1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6 7 7 7 8 8 8 9 9 9 10 10 10];
% xx=[1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2 2 1 1 1 2 2];
% ob(45:85,25:65)=repmat(xx,[41,1]);
% ob=repmat(1:32,[3,1]);
% ob=ob(:);
% ob=repmat(ob',[108,1]);


% ob(45:75,25:55)=1;
% ob(31:94,21:84)=phantom(64);
% ob=[1:96];
% ob=repmat(ob,[140 1]);

% ob=phantom(64);
% ob(18:45,18:45)=ob(18:45,18:45)+0.5;
% ob=padarray(ob,[16 16 0],'both');
% obj = bsxfun(@times, ob, sens);
% obj=ob;

% obj = bsxfun(@times, sum(obj.*conj(sens),3), sens);


rng(10000)
t=1:1024;
% off_x=[0 pi/1 pi/2 pi/3 pi/4 pi/5];
off_x=0;
slo=randn(1,5)/2;
off_y=randn(1,5)/2;
% for i=1:5
%     ww(i,:)=sin(t*(pi/30)+off_x)*slo(i)+off_y(i);
% end

ww=ww(:,1:frame);
% idx=randperm(frame,frame);
% ww_new=ww(:,idx);
% ww=ww_new;

% ww=ww(:,1:2:end);

% Set Magnitude
% obj = phantom(64);
% par.kx=64;
% par.ky=64;
% par.shot=4;
% obj = bsxfun(@times, obj, sens);
% obj=obj(7:102,:,:);
[par.kx,par.ky,par.coil]=size(obj);
par.kx=128;par.ky=128;

% Set 4 different smooth phase profiles
% [x, y] = meshgrid(0:63, 0:63);
% phs(:,:,1) = exp(1j*2*pi*(x)/64);
% phs(:,:,2) = exp(1j*2*pi*(y)/32);
% phs(:,:,3) = exp(1j*2*pi*(x+2*y)/48);
% phs(:,:,4) = exp(1j*2*pi*(3*x-y)/16);
% phs(:,:,2) = exp(1j*2*pi*(x)/30/10);
% phs(:,:,3) = exp(1j*2*pi*(x)/45/10);
% phs(:,:,4) = exp(1j*2*pi*(x)/20)/10;

np=5;
% par.kx=128
% par.ky=128
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

[~,frame]=size(ww);
phs=zeros(par.kx,par.ky,frame);
for i=1:frame
for j=1:np
    phs(:,:,i)=phs(:,:,i)+p(:,:,j)*ww(j,i);
end
end
phs=permute(repmat(phs,[1,1,1,par.coil]),[1,2,4,3]);

% obj= sum(obj.*conj(par.sens),3);%n*l*coil*shot
% obj=abs(obj);
% load('smoothp.mat')
% obj=obj.*exp(1j*p);
% obj = bsxfun(@times, obj, par.sens);


% Generate complex obj
obj = bsxfun(@times, obj, exp(1j*phs*1));
% obj=repmat(obj,[1 1 1 frame]);
 
 
% %  [ traj,par.dcf ] = calc_traj( ones(108, 96) ,'sequential',0,[0.2 -0.3]);
%  [ traj,par.dcf ] = calc_traj( ones(140, 96),'golden',0,[0 0]);
%  %% 
% par.FT = NUFFT(traj,[],1,[0 0],[140 96],2);
% % for i=1:frame
% %     for j=1:par.coil
% %         k_obj(:,:,j,i)=par.FT*obj(:,:,j,i);
% %     end
% % end
%  par.MCFT = MCNUFFT(traj,[],sens);
%   
% for i=1:frame
% %     for j=1:par.coil
%         k_obj(:,:,:,i)=par.MCFT*sum(obj(:,:,:,i).*conj(sens),3);
% %     end
% end


% generate k-space
k_obj   = fftdim(obj,1:2);
% par.k_gt=k_obj;

% rng(t+100);
% noise_i = randn(size(k_obj))*noise_level;
% noise_j = randn(size(k_obj))*noise_level;
% k_obj=noise_i+1j*noise_j+k_obj;

%  mest=sum(ww,1).';
%  mest=repmat(mest,[1,140,140]);
%  mest=permute(mest,[2 3 1]);
% 
%  
% kgrid=repmat((-70:69)'/140,[1,96]);
% kgrid=repmat(kgrid,[1 1 2048]);
% k_obj=k_obj.*exp(1j*kgrid.*mest);
% %   m=bsxfun(@times,repmat(kgrid,[1 1 2048]),mest);


end