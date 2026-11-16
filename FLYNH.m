function [f_result,rr]=FLYNH(p,c,N)
% fit
cn=[0.5 1:N];%wave number
sigma=2*pi/(p(end)-p(1));
n=length(cn);
for j=1:n
    fu(j)=cn(j)*sigma;
end

tc=zeros(length(p),2*n+1,'double');
tc(:,1)=1;
for w=1:length(p)
    tc(w,2:n+1)=cos(p(w)*fu);
    tc(w,n+2:2*n+1)=sin(p(w)*fu);
end
[A,~,rr,~,~]=regress(c,tc);
f_result=tc*A;

end