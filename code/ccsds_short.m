function H = ccsds_short(len)
% CCDS codes taken from the orange book "SHORT BLOCK LENGTH LDPC CODES FOR TC
% SYNCHRONIZATION AND CHANNEL CODING", April 2015: CCSDS 231.1-O-1

switch len
    case 128
        Hc = [-1 2 14 6 -1 0 13 0
            6 -1 0 1 0 -1 0 7
            4 1 -1 14 11 0 -1 3
            0 1 9 -1 14 1 0 -1];
        c = 128/8;
        H = ldpc.utils.gen_QC_H(Hc,c);
        H(1:c,1:c) = mod(speye(c) + circshift(speye(c),[0 7]),2);
        H(c+1:2*c,c+1:2*c) = mod(speye(c) + circshift(speye(c),[0 15]),2);
        H(2*c+1:3*c,2*c+1:3*c) = mod(speye(c) + circshift(speye(c),[0 15]),2);
        H(3*c+1:4*c,3*c+1:4*c) = mod(speye(c) + circshift(speye(c),[0 13]),2);
    case 256
        Hc = [-1 15 25 0 -1 20 12 0
            28 -1 29 24 0 -1 1 20
            8 0 -1 1 29 0 -1 21
            18 30 1 -1 25 26 0 -1];
        c = 256/8;
        H = ldpc.utils.gen_QC_H(Hc,c);
        H(1:c,1:c) = mod(speye(c) + circshift(speye(c),[0 31]),2);
        H(c+1:2*c,c+1:2*c) = mod(speye(c) + circshift(speye(c),[0 30]),2);
        H(2*c+1:3*c,2*c+1:3*c) = mod(speye(c) + circshift(speye(c),[0 28]),2);
        H(3*c+1:4*c,3*c+1:4*c) = mod(speye(c) + circshift(speye(c),[0 30]),2);
    case 512
        Hc = [-1 30 50 25 -1 43 62 0
            56 -1 50 23 0 -1 37 26
            16 0 -1 27 56 0 -1 43
            35 56 62 -1 58 3 0 -1];
        c = 512/8;
        H = ldpc.utils.gen_QC_H(Hc,c);
        H(1:c,1:c) = mod(speye(c) + circshift(speye(c),[0 63]),2);
        H(c+1:2*c,c+1:2*c) = mod(speye(c) + circshift(speye(c),[0 61]),2);
        H(2*c+1:3*c,2*c+1:3*c) = mod(speye(c) + circshift(speye(c),[0 55]),2);
        H(3*c+1:4*c,3*c+1:4*c) = mod(speye(c) + circshift(speye(c),[0 11]),2);
    otherwise
        error('Blocklength %d is not defined in the standard.\n', len); 
end




end