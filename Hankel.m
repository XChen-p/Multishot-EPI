function [H_ms] =Hankel(k_ms,r,dimension)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

 % r: filter size

H_ms=[];

if dimension =='3D'
[m,n,l,coil,shot]=size(k_ms);

for s=1:shot
    k=k_ms(:,:,:,:,s);
    H_ss = zeros((m-r+1)*(n-r+1),r^2,coil);
    for c=1:coil
        for q=1:l-r+1
            for j=1:n-r+1
                for i=1:m-r+1
                temp=k(i:i+r-1,j:j+r-1,q:q+r-1,c);
                inx=i+(j-1)*(m-r+1)+(q-1)*(m-r+1)*(n-r+1);
                H_ss(inx,:,c)=temp(:).';
                end
            end
        end
    end  
    H_ms=[H_ms reshape(H_ss,[],r*r*r*coil)]; 
end

else
    if ndims(k_ms)==3
        k_ms = permute(k_ms,[1,2,4,3]);% m*n*coil*shot
    end
    [m,n,coil,shot]=size(k_ms);
for s=1:shot
    k=k_ms(:,:,:,s);
    H_ss = zeros((m-r+1)*(n-r+1),r^2,coil);
    for c=1:coil
        for j=1:n-r+1
            for i=1:m-r+1
                temp=k(i:i+r-1,j:j+r-1,c);
                inx=i+(j-1)*(m-r+1);
                H_ss(inx,:,c)=temp(:).';
            end
        end
    end
     
    H_ms=[H_ms reshape(H_ss,[],r*r*coil)]; 
end
    
 
end

end

