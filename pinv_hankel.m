function [M,N] = pinv_hankel(H, k, Nx, Ny, Nc,Nshot)
    M = zeros(Nx*Ny, Nc,Nshot);
    N = zeros(Nx*Ny, Nc,Nshot);
    H = reshape(H,Nx-(k-1),[],k,k,Nc,Nshot);
    for ii = 1:size(H,2)
        for jj = 1:k
            for kk = 1:k
                idx = ((ii-1)*Nx+1:(ii-1)*Nx+Nx-(k-1))+(jj-1)+(kk-1)*Nx;
                M(idx,:,:) = M(idx,:,:) + permute(H(:,ii,jj,kk,:,:),[1,5,6,2,3,4]);
                N(idx,:,:) = N(idx,:,:) + 1;
            end
        end
    end
    M = reshape(M./N, Nx,Ny,Nc,Nshot);    
    N = reshape(N, Nx,Ny,Nc,Nshot); 
end
