function  H = gen_QC_H(Hc,c)
% generate quasi-cyclic parity check matrix of an array of size(Hc) cyclic
% shifted identity matrices of size c

ID = speye(c);
[e,f] = size(Hc);

[row_idx,col_idx] = find(Hc~=-1);
nnz = length(row_idx)*c*e*f;

H = spalloc(e*c,f*c,nnz);

for ii=1:length(row_idx)
    H((row_idx(ii)-1)*c+1:row_idx(ii)*c,(col_idx(ii)-1)*c+1:col_idx(ii)*c) = circshift(ID,[0 Hc(row_idx(ii),col_idx(ii))]);    
end

end