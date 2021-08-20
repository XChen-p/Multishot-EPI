
% Main Script
filepath ='/Users/xchen/Documents/MATLAB/Ox/20200224/d_kobj_18mm_8_2_Rz1/';
file = dir([filepath 'd_kobj_18mm_8_2_Rz1*.mat']);

load('sens_d_18mm_espiritNew_r96t3c6.mat')
par.sens=sens;

acc=2;%prospective/retrospective R_3D
for t=1:length(file)
    load([filepath file(t).name])
    
    [par.kx,par.ky,par.coil,par.shot]=size(k_obj);
     par.shot=par.ky/acc;% retrospective undersampling
    
    % sampling trajectory
    sample = squeeze(CAIPI_Sampling(1,[1,par.kx,par.ky,par.shot], 8,2,2, 1,0));
    par.sample = permute(repmat(sample,[1,1,1,par.coil]),[1,2,4,3]);
    
    
    % retrospective undersampling
    k_obj=k_obj(:,:,:,1:par.shot);
    k_obj = k_obj.*par.sample;
    
    % binning into segments
    par.state=8/acc;
    [k_obj,par] = shot2state(k_obj,par.state,par,'sequential');
    
    % shot-combined SENSE recon
    output1 = MCsense(sum(k_obj(:,:,:,1:end),4),sum(par.sample(:,:,:,1:end),4),par.sens);
    
    figure(1)
    imshow(abs(output1),[])
    img1(:,:,t)=output1;

    % structured low-rank recon
    par.coilcom=1;


    par.f       = 6;
    par.lambda  = 1E-4;
    par.rho     = par.lambda*1E-2;
    par.niter   = 10;
    % rho adjustment parameter
    par.m = 2;
    % over-relaxation parameter
    par.r = 1.5;

    % turn on or off printing cost function
    % printing cost function can slow things down
    par.verbose = 1;

    % Choose initialisation for x

    x0 = repmat(sum(k_obj,4)./(sum(par.sample,4)+eps),[1,1,1,par.shot]);
    x0    = M_adj(x0,par);
    % solve ADMM
    output2 = ADMM(k_obj, x0, par);
    figure(2)
    imshow(abs(output2),[])
    img2(:,:,t)=output2;


end


function k_sc = M_adj(k_mc, par)
im_mc = ifftdim(k_mc,1:2);%n*l*coil*shot
im_sc = sum(im_mc.*conj(par.sens),3);%n*l*1*shot
k_sc  = fftdim(im_sc,1:2);
end


