function [k,par] = shot2state(k,state,par,mode)


[kx,ky,coil,shot]=size(k);
num=shot/state;
k_mst=zeros(kx,ky,coil,state);
sample_mst=zeros(kx,ky,coil,state);


if strcmp(mode,'sequential')
    for s=1:state
        k_mst(:,:,:,s) = sum(k(:,:,:,1+(s-1)*num:s*num),4);
        sample_mst(:,:,:,s) = sum(par.sample(:,:,:,1+(s-1)*num:s*num),4);
    end
else strcmp(mode,'respiratory') % not applicable
    
    for s=1:state
        idx=sort(par.idx(1+(s-1)*num:s*num));
        sample = par.sample(:,:,:,idx);
        kk = k(:,:,:,idx);
        k_mst(:,:,:,s) = sum(kk,4);
        sample_mst(:,:,:,s) = sum(sample,4);
 
    end
         
end
 

k = k_mst;
par.sample = sample_mst;
par.shot=state;
end

