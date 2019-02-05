function info = get_H_info(H)

info.nc = size(H,2);
info.mc = size(H,1);

% we assume this for the moment
info.nct = info.nc;
info.mct = info.mc;
info.puncture = [];
info.shorten = [];

info.nnz = full(sum(H(:)));
[tmp1,tmp2] = find(H');
info.r = int32(tmp2-1);
info.c = int32(tmp1-1);

info.rw = full(sum(H,2));
info.rw = int32(info.rw(:));
info.cw = full(sum(H,1));
info.cw = int32(info.cw(:));

MAX_ROW_W = 40;
MAX_COL_W = 40;

% z = 0;
% for ii=1:info.M
%     idx = find(H(ii,:));
%     info.rn((ii-1)*MAX_ROW_W+1:(ii-1)*MAX_ROW_W+length(idx)) = z:(z+length(idx)-1);
%     z = z + length(idx);
% end
info.rn = zeros(MAX_ROW_W*info.mc,1);
w = zeros(info.mc,1);
for ii=1:info.nnz
    info.rn(info.r(ii)*MAX_ROW_W+w(info.r(ii)+1)+1) = ii-1;
    w(info.r(ii)+1) = w(info.r(ii)+1)+1;
end
info.rn = int32(info.rn(:));

% z = 0;
% for ii=1:info.N
%     idx = find(H(:,ii));
%     info.cn((ii-1)*MAX_ROW_W+1:(ii-1)*MAX_COL_W+length(idx)) = z:(z+length(idx)-1);
%     z = z + length(idx);
% end

% info.cn = zeros(MAX_COL_W*info.N);
% for ii=1:info.N
%     for jj=1:info.cw(ii)
%         [~,info.cn((ii-1)*MAX_COL_W+jj)] = ismember([jj-1 ii-1],[info.r info.c], 'rows');
%     end
% end
info.cn = zeros(MAX_COL_W*info.nc,1);
w = zeros(info.nc,1);
for ii=1:info.nnz
    info.cn(info.c(ii)*MAX_COL_W+w(info.c(ii)+1)+1) = ii-1;
    w(info.c(ii)+1) = w(info.c(ii)+1)+1;
end
info.cn = int32(info.cn(:));

end
