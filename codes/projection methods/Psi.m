function C=Psi(b,p_c)

order=length(p_c)-1;
C=p_c*chebpol(order,b)'./(sum (abs(p_c)));

end