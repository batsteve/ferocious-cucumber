function dy=oderhs(~,y, par)
global n1 n2

% t
u1=reshape(y(1:n1*n2),        n2,n1);
u2=reshape(y(n1*n2+1:2*n1*n2),n2,n1);
[rhs1, rhs2]=rhs(u1,u2, par);
rhs1_vec = reshape(rhs1,n1*n2,1);
rhs2_vec = reshape(rhs2,n1*n2,1);
dy = [rhs1_vec;rhs2_vec];

