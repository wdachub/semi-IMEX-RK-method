function D=diffmat_period(x,T,n_s,m)
%x grid points
%T period
%n_s number of stencil points
%m required order of derivatives
% N=length(x);
n2=ceil(n_s/2);
x2=[x(end-n2+1:end)-T+x(1); x; T+x(1:n2)];
D_aux=diffmat(x2,n_s,m);

% D_aux{1}=full(D_aux{1});
for i=m:-1:1
D{i}=D_aux{i}(1+n2:end-n2,1+n2:end-n2);
D{i}(1:n2,end-n2+1:end)=D{i}(1:n2,end-n2+1:end)+D_aux{i}(n2+1:2*n2,1:n2);
D{i}(end-n2+1:end,1:n2)=D{i}(end-n2+1:end,1:n2)+D_aux{i}(end-2*n2+1:end-n2,end-n2+1:end);


% D{i}(1:n2-1,end-n2+2:end)=D{i}(1:n2-1,end-n2+2:end)+D_aux{i}(n2+1:2*n2-1,2:n2);
% D{i}(end-n2+2:end,1:n2-1)=D{i}(end-n2+2:end,1:n2-1)+D_aux{i}(end-2*n2+2:end-n2,end-n2+1:end-1);

%%
% 
%  PREFORMATTED
%  TEXT
% 

% 1
end


   