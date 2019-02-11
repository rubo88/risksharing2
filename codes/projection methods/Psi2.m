function C=Psi2(b,p_b)

order=length(p_b)-1;
C=p_b*chebpol(order,b)'./(sum (abs(p_b)));

end