function [Wk] = removeSmallValues(W_k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

for q = 1:length(W_k)
    if abs(real(W_k(q))) < 1e-3
        W_Re(q) = 0;
    else W_Re(q) = real(W_k(q));
    end;
    
    if abs(imag(W_k(q))) < 1e-3
        W_Im(q)=imag(W_k(q))*0;
    else W_Im(q)=imag(W_k(q))
    end;
    Wk(q)=W_Re(q)+i*W_Im(q);
end;


end

