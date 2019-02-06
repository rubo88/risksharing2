function q=Phi(b,p_q)
order=length(p_q)-1;
q=p_q*chebpol(order,b)';
end