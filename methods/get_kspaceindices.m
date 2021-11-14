%%% define k-space indices
function [ind_full,ind_filter,ind_filter2,k] = get_kspaceindices(overres,res,filter_siz)

k = gpuArray.zeros(3,prod(overres));
n1 = overres(1);n2 = overres(2);n3 = overres(3);
v1_o = mod(n1,2);
v2_o = mod(n2,2);
v3_o = mod(n3,2);

if v2_o==0   %%% even
    kx = ifftshift(repmat(-n2/2:(n2/2)-1,n1,1));
else
    kx = ifftshift(repmat(-(n2-1)/2:(n2-1)/2,n1,1));
end

if v1_o==0  %%% even
    ky = ifftshift(repmat((-n1/2:(n1/2)-1).',1,n2));
else
    ky = ifftshift(repmat((-(n1-1)/2:(n1-1)/2).',1,n2));
end

kx = repmat(kx,1,1,n3);
ky = repmat(ky,1,1,n3);

if v3_o==0 %%% even
    temp = [0:(n3/2)-1 , -n3/2 :-1];
else
    temp = [0:(n3-1)/2 , -(n3-1)/2 :-1];
end

for i=1:n3
    kt(:,:,i) = temp(i)*gpuArray.ones(n1,n2);
end

k(1,:) = kx(:);
k(2,:) = ky(:);
k(3,:) = kt(:);

% clear kx ky kt

na = res(1);nb = res(2); nc = res(3);
v1_r = mod(na,2);
v2_r = mod(nb,2);
v3_r = mod(nc,2);

if v3_r==0
    if (v2_r == 0) && (v1_r == 0)
        ind_full = find(((k(1,:) >= -nb/2) & (k(1,:) < nb/2)) & ((k(2,:) >= -na/2) & (k(2,:) < na/2)) & ((k(3,:) >= -nc/2) & (k(3,:) < nc/2)));
    elseif (v2_r ==1) && (v1_r == 0)
        ind_full = find((abs(k(1,:)) <= (nb-1)/2) & ((k(2,:) >= -na/2) & (k(2,:) < na/2)) & ((k(3,:) >= -nc/2) & (k(3,:) < nc/2)));
    elseif (v2_r==0) && (v1_r==1)
        ind_full = find(((k(1,:) >= -nb/2) & (k(1,:) < nb/2)) & (abs(k(2,:)) <= (na-1)/2) & ((k(3,:) >= -nc/2) & (k(3,:) < nc/2)));
    else
        ind_full = find((abs(k(1,:)) <= (nb-1)/2) & (abs(k(2,:)) <= (na-1)/2) & ((k(3,:) >= -nc/2) & (k(3,:) < nc/2)));
    end
else
    if (v2_r == 0) && (v1_r == 0)
        ind_full = find(((k(1,:) >= -nb/2) & (k(1,:) < nb/2)) & ((k(2,:) >= -na/2) & (k(2,:) < na/2)) & (abs(k(3,:)) <= (nc-1)/2));
    elseif (v2_r ==1) && (v1_r == 0)
        ind_full = find((abs(k(1,:)) <= (nb-1)/2) & ((k(2,:) >= -na/2) & (k(2,:) < na/2)) & (abs(k(3,:)) <= (nc-1)/2));
    elseif (v2_r==0) && (v1_r==1)
        ind_full = find(((k(1,:) >= -nb/2) & (k(1,:) < nb/2)) & (abs(k(2,:)) <= (na-1)/2) & (abs(k(3,:)) <= (nc-1)/2));
    else
        ind_full = find((abs(k(1,:)) <= (nb-1)/2) & (abs(k(2,:)) <= (na-1)/2) & (abs(k(3,:)) <= (nc-1)/2));
    end
end

fa = filter_siz(1); fb = filter_siz(2); fc = filter_siz(3);
v1_f = mod(fa,2);
v2_f = mod(fb,2);
v3_f = mod(fc,2);

if v3_f == 0
    
    if (v2_f == 0) && (v1_f == 0)
        ind_filter = find(((k(1,:) >= -fb/2) & (k(1,:) < fb/2)) & ((k(2,:) >= -fa/2) & (k(2,:) < fa/2)) & ((k(3,:) >= -fc/2) & (k(3,:) < fc/2)));
    elseif (v2_f ==1) && (v1_f == 0)
        ind_filter = find((abs(k(1,:)) <= (fb-1)/2) & ((k(2,:) >= -fa/2) & (k(2,:) < fa/2)) & ((k(3,:) >= -fc/2) & (k(3,:) < fc/2)));     
    elseif (v2_f==0) && (v1_f==1)
        ind_filter = find(((k(1,:) >= -fb/2) & (k(1,:) < fb/2)) & (abs(k(2,:)) <= (fa-1)/2) & ((k(3,:) >= -fc/2) & (k(3,:) < fc/2)));
    else
        ind_filter = find((abs(k(1,:)) <= (fb-1)/2) & (abs(k(2,:)) <= (fa-1)/2) & ((k(3,:) >= -fc/2) & (k(3,:) < fc/2)));
    end
    
else
    
    if (v2_f == 0) && (v1_f == 0)
        ind_filter = find(((k(1,:) >= -fb/2) & (k(1,:) < fb/2)) & ((k(2,:) >= -fa/2) & (k(2,:) < fa/2)) & (abs(k(3,:)) <= (fc-1)/2));
    elseif (v2_f ==1) && (v1_f == 0)
        ind_filter = find((abs(k(1,:)) <= (fb-1)/2) & ((k(2,:) >= -fa/2) & (k(2,:) < fa/2)) & (abs(k(3,:)) <= (fc-1)/2));     
    elseif (v2_f==0) && (v1_f==1)
        ind_filter = find(((k(1,:) >= -fb/2) & (k(1,:) < fb/2)) & (abs(k(2,:)) <= (fa-1)/2) & (abs(k(3,:)) <= (fc-1)/2) );
    else
        ind_filter = find((abs(k(1,:)) <= (fb-1)/2) & (abs(k(2,:)) <= (fa-1)/2) & (abs(k(3,:)) <= (fc-1)/2));
    end
    
end
    
ind_filter2 = find((abs(k(1,:)) <= (fb-1)) & (abs(k(2,:)) <= (fa-1)) & (abs(k(3,:)) <= (fc-1)));

end
        
        
    
    

