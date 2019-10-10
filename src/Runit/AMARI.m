function a = AMARI( W,CC,P)% function a = AMARItest( G)%% Returns a performance measure based on the system’s global matrix using the AMARI criterion% Cichocki, A and Amari, S "Adaptive Blind Signal and Image Processing"% John Wiley & Sons 2005 p308
G=W*CC*P;
a = sum(sum(abs(G),2)./max(abs(G),[],2)-1) + sum(sum(abs(G))./max(abs(G))-1);
a = a/size(G,1);