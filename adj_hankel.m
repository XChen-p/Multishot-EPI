function [M] = adj_hankel(H, r, Nx, Ny, Nc,Nshot)
num=(Nx-r+1)*(Ny-r+1);
M=zeros(Nx,Ny,Nc,Nshot);
for i=1:num
    idm=mod(i-1,(Nx-r+1))+1;
    row=idm:idm+r-1;
    idn=floor((i-1)/(Nx-r+1))+1;
    column=idn:idn+r-1;

    M(row,column,:,:)=M(row,column,:,:)+reshape(H(i,:),[r r Nc Nshot]);
end

end