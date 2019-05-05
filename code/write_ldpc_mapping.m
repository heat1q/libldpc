function write_ldpc_mapping(filename,Hinfo,mode)

% shortened and punctured bits are not transmitted
% but punctured bits may be attributed to one bit level anyway
bits_pos = 0:Hinfo.nc-1; % watch out for proper indexing (zero start ...)
bits_pos = setdiff(bits_pos, Hinfo.shorten);
bits_pos = setdiff(bits_pos, Hinfo.puncture); % watch out for proper indexing (zero start ...)

mapping = zeros(mode.bits, mode.n);

if strcmp(mode.name, 'random')
    idx = randperm(length(bits_pos));
    idx = idx(:);
    idx = reshape(idx,mode.bits,[]);
    for ii=1:mode.bits
        mapping(ii,:) = bits_pos(idx(ii,:));
    end  
elseif strcmp(mode.name, 'ordered')
    if ~isfield(mode,'order')
        error('missing parameter "order"');
    end
    if isfield(mode, 'order') && length(mode.order) ~= mode.bits
        error('length of order vector does not fit')
    end
    if isfield(mode, 'order') && ~isequal(unique(mode.order), 1:mode.bits)
        error('order vector does not make sense')
    end
    idx = bits_pos(:);
    idx = reshape(idx, [],  mode.bits);
    for ii=1:mode.bits
        mapping(mode.order(ii),:) = idx(:, ii);
    end
elseif strcmp(mode.name, 'consecutive')
    idx = bits_pos;
    idx = idx(:);
    idx = reshape(idx, mode.bits, []);
    for ii=1:mode.bits
        mapping(ii,:) = idx(ii,:);
    end
elseif strcmp(mode.name, 'manual')
    mapping = mode.mapping;
else
    error('desired mode not known');
end

f = fopen(filename, 'w');
for ii=1:mode.bits
    fprintf(f, '%d, ', mapping(ii,1:end-1));
    fprintf(f, '%d\n', mapping(ii,end));
end
fclose(f);

end