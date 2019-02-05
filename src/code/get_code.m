n = 2400;
m = 600;
nnz = 7200;

table = readtable('code_3x12.txt');

cw = sparse([],[],[],1,n,0);
index = [1245 1540 1602 1760 2067 2083 2261 2325];


for ii=1:length(index)
    cw(index(ii)+1) = 1;
end

m_j = table.J;
n_i = table.I;

H = sparse([],[],[],m,n,0);

for ii=1:nnz
	H(m_j(ii)+1,n_i(ii)+1) = 1;
end

check = sum(cw*H')