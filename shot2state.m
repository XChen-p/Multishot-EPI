function [k,par] = shot2state(k,state,par,mode)


[kx,ky,coil,shot]=size(k);
num=shot/state;
k_mst=zeros(kx,ky,coil,state);
sample_mst=zeros(kx,ky,coil,state);

% x=-5.5:5.5;
% % x=-11:0;
% y=normpdf(x,0,8);
% y=y./max(y);
% figure
% plot(y)
% y=sqrt(y);
% figure
% plot(y)
% y2=y(end:-1:1);
% w=ones(kx,ky,coil,num);
% for i=1:12
% w1(:,:,:,i)=w(:,:,:,i)*y(i);
% w2(:,:,:,i)=w(:,:,:,i)*y2(i);
% end


if strcmp(mode,'sequential')
    for s=1:state
%         if mod(s,2)==0
%             w=w2;
%         else
%             w=w1;
%         end
        k_mst(:,:,:,s) = sum(k(:,:,:,1+(s-1)*num:s*num),4);
        sample_mst(:,:,:,s) = sum(par.sample(:,:,:,1+(s-1)*num:s*num),4);
 
%             figure(10)
%             imshow(squeeze(sample_mst(:,:,1,s)))
    end
else strcmp(mode,'respiratory')
    
    for s=1:state
        idx=sort(par.idx(1+(s-1)*num:s*num));
        sample = par.sample(:,:,:,idx);
        kk = k(:,:,:,idx);
        k_mst(:,:,:,s) = sum(kk,4);
        sample_mst(:,:,:,s) = sum(sample,4);
                    figure(10)
            imshow(squeeze(sample_mst(:,:,1,s)))
    end
        
        
 
%     k = k(:,:,:,par.idx);
%     par.sample = par.sample(:,:,:,par.idx);
%     for s=1:state
%         k_mst(:,:,:,s) = sum(k(:,:,:,1+(s-1)*num:s*num),4);
%         sample_mst(:,:,:,s) = sum(par.sample(:,:,:,1+(s-1)*num:s*num),4);
%         
%             figure(10)
%             imshow(squeeze(sample_mst(:,:,1,s)))
%     end 
    
end
figure(20)
sam=squeeze(sum(sample_mst(:,:,1,:),4));
imshow(sam,[])


k = k_mst;
par.sample = sample_mst;
par.shot=state;
end

