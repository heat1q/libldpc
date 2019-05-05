function write_QC_layerfile(filename, B, L)
%Write layers of QC code constructed from base matrix B and factor L

f = fopen(filename,'w');

sizeB = size(B); 
fprintf(f, 'nl: %d\n', sizeB(1));
for ii=0:sizeB-1
    fprintf(f, 'cn[i]: %d\n', L);
    for jj=(ii*L):(ii+1)*L-1
        fprintf(f, '%d\n', jj);
    end
end
end

