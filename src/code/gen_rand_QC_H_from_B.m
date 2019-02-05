function H = gen_rand_QC_H_from_B(B,L)
% generate QC 
% with random cyclic shift, based on entry in B
[M,N] = size(B);

H = sparse(M*L,N*L);

for ii=1:M
    for jj=1:N
        if B(ii,jj) ~= 0
            idx = randperm(B(ii,jj),1);
            
            tmp = circshift(speye(L),[1, idx]);
            
            H((ii-1)*L+1:ii*L,(jj-1)*L+1:jj*L) = tmp;
        end
    end
end
end

