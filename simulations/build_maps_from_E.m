function [T1, T2, PD, df] = build_maps_from_E(E, LUT, N1, N2)
if nargin < 3
    % Assumption: E has N*N lines
    N1 = sqrt(size(E, 1));
    N2 = N1;
end

T1 = zeros(N1, N2, 1, 'single');
T2 = T1;
df = T1;
PD = single(full(abs(sum(E, 2))));

ind = find(PD);
% for i=1:length(ind)
%     res=(abs(E(ind(i),:))*LUT)/sum(abs(E(ind(i),:)));
%     T1(ind(i))=res(1);
%     T2(ind(i))=res(2);
%     if (size(LUT,2)==3)
%         df(ind(i))=res(3);
%     end
% end

res = (E > 0) * LUT;
T1(ind) = res(ind, 1);
T2(ind) = res(ind, 2);
if (size(LUT, 2) == 3)
    df(ind) = res(ind, 3);
end

% T1=reshape(T1,N1,N2);
% T2=reshape(T2,N1,N2);
PD = reshape(PD, N1, N2);
% df=reshape(df,N1,N2);

end
