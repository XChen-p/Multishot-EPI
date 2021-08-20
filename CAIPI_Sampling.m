function mask = CAIPI_Sampling(verbose,dims, width, slope, R, pf_y,pf_z)
%
%   dims is [Nx, Ny, Nz, Nt]
%   width is the extent of the "band" in kz for blipping
%   slope is the dkz/dky for blipping
%   shift is the "block" offset for each shot
%   R is the in-plane acceleration factor
%   pf_y is the partial fourier factor in ky (1 for full sampling)
%   pf_z is the partial fourier factor in kz (1 for full sampling)
%


%   check if dims(3)/width is an even number
if mod(dims(3),width)
    error('dims(3) must be divisible by width');
end

if pf_z~=1
    dims(4)=ceil(dims(4)*pf_z);
end


%   generate "shift" parameters

num_block = ceil(dims(3)/width * pf_z);

%   shift_1 is the "between block" shift

% %   use golden ratio, but if not a good value, increment until good
% shift_1 =   round((num_block)*(sqrt(5)-1)/2);
% while lcm(shift_1, num_block) ~= shift_1*num_block
%     shift_1 = shift_1 + 1;
% end
%use sequential order
shift_1=1;

%   shift_2 is the "within block" shift
%golden ratio
shift_2 =   round(width*(sqrt(5)-1)/2);
while lcm(shift_2, width) ~= shift_2*width
    shift_2 = shift_2 + 1;
end
% shift_2=dims(3)/dims(4);




%   initialise sampling mask
% mask    =   false(dims);
mask    =   zeros(dims);

%   loop over time points (shots) and ky
% figure
for t = 1:dims(4)
    
    idx_s =   floor((t-1)/(num_block));
    
    y0  =   mod(t-1, R);
    %     y0  =   mod(mod(t-1,num_block), R);%no difference if mod(num_block,R)=0
    z0  =   mod(shift_1*(t-1), num_block)*width;
    b0  =   mod(shift_2*idx_s, width);
    
    for i = ceil(dims(2)*(1-pf_y))+1:R:dims(2)
        %         ii = (ceil(i/R)-1)*slope+1;%
        ii=(i-1)*slope+1;
        mask(:, y0 + i, z0 + mod(b0 + ii - 1, width) + 1, t) =  true;
    end
    if verbose==1
        figure(2)
        s=sum(squeeze(mask(1,:,:,:)),3);
        imshow(s,[],'border','tight','colormap',jet)
    end
    
end
