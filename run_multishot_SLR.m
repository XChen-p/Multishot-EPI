
% Main Script

 
filepath ='/data/';
% file = dir([filepath 'slice64_2_1_*.mat']); % a 2D slice from the 3D dataset shown in Fig. 4, seg-CAIPI(2,1), R=2x2
file = dir([filepath 'slice64_8_3_*.mat']); % a 2D slice from the 3D dataset shown in Fig. 4, seg-CAIPI(8,3), R=2x2

load('data/slice64_sens.mat')
% load('sens_1_5.mat')% coil sensitivity for simulation data 
par.sens=sens;

acc=2;%prospective/retrospective R_3D
for t=1:length(file)
    load([filepath file(t).name])% in vivo   
       %     [k_obj] = gen_obj(t,48,3E-6); % generate the simulation data
    [par.kx,par.ky,par.coil,par.shot]=size(k_obj);
    
    % sampling trajectory
    sample = squeeze(CAIPI_Sampling(0,[1,par.kx,par.ky,par.shot/acc],8,3,2,1,1));
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
    par.lambda  = 3E-4;
    par.rho     = par.lambda*1E-2;
    par.niter   = 10;
    % rho adjustment parameter
    par.m = 2;
    % over-relaxation parameter
    par.r = 1;%1.5;

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


