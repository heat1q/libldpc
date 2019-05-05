function write_ldpc_codefile(filename, info, do_encode)

f = fopen(filename,'w');

fprintf(f, 'nc: %d\n', info.nc);
fprintf(f, 'mc: %d\n', info.mc);
fprintf(f, 'nct: %d\n', info.nct);
fprintf(f, 'mct: %d\n', info.mct);
fprintf(f, 'nnz: %d\n', info.nnz);
if ~isempty(info.puncture)
    fprintf(f, 'puncture [%d]: ', length(info.puncture));
    fprintf(f, '%d ', info.puncture(1:end-1));
    fprintf(f, '%d\n', info.puncture(end));
else
    fprintf(f, 'puncture [0]: \n');
end
if ~isempty(info.shorten)
    fprintf(f, 'shorten [%d]: ', length(info.shorten));
    fprintf(f, '%d ', info.shorten(1:end-1));
    fprintf(f, '%d\n', info.shorten(end));
else
    fprintf(f, 'shorten [0]: \n');
end

for ii=1:info.nnz
    fprintf(f,'%d %d\n', info.r(ii), info.c(ii));
end

if do_encode
    H = sparse(double(info.r)+1,double(info.c)+1,1);
    G = ldpc.utils.get_G(H);
    [k ,n] = size(G);
    
    genmat = cell(n-k,1);
    for ii=k+1:n
        genmat{ii-k} = [sum(G(:,ii))  find(G(:,ii)).'-1];
    end
    
    for ii=1:(n-k)
        fprintf(f, '%d ', genmat{ii});
        fprintf(f, '\n');
    end
end

fclose(f);

end
