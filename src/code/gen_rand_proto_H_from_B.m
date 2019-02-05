function H = gen_rand_proto_H_from_B(B,L)
% generate protograph
% with random shift & add, based on entry in B

[M,N] = size(B);

H = sparse(M*L,N*L);

for ii=1:M
    for jj=1:N
        if B(ii,jj) ~= 0
            idx = randperm(L,B(ii,jj));
            tmp = 0;
            for kk=1:B(ii,jj)
                tmp = tmp + circshift(speye(L),[1, idx(kk)]);
            end
            
            H((ii-1)*L+1:ii*L,(jj-1)*L+1:jj*L) = tmp;
        end
    end
end

end