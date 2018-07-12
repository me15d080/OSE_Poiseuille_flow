% Read data:
file1  = '/home/vikas/Desktop/hw_5430/Ar';
file2  = '/home/vikas/Desktop/hw_5430/Ai';
file3  = '/home/vikas/Desktop/hw_5430/Br';
file4  = '/home/vikas/Desktop/hw_5430/Bi';


Ar = load(file1,'-ascii');
Ai = load(file2,'-ascii');
Br = load(file3,'-ascii');
Bi = load(file4,'-ascii');

N  = sqrt(length(Ar));

AR =zeros(N,N);
AI =zeros(N,N);
BR =zeros(N,N);
BI =zeros(N,N); 
for i=1:N 
    for j=1:N
        AR(i,j)=Ar((i-1)*N+j);
        AI(i,j)=Ai((i-1)*N+j);
        BR(i,j)=Br((i-1)*N+j);
        BI(i,j)=Bi((i-1)*N+j); 
    end 
end
i=sqrt(-1);
A = AR+i*AI
B = BR+i*BI

% find eigenvalues 
[V,D]=eig(A,B);
for j=1:N;
real_c(j)=real(D(j,j));
imag_c(j)=imag(D(j,j));
end
%--Eigen Modes--%
psi_hat=zeros(1,N);
y = linspace(-1,1,N);
for m=N-2:-1:2;
yC=cos(m*pi/(N-1));
cheb=zeros(1,N);% contains Chebyschev polynomial and derivs
% Cheb Poly @ yC
cheb(1,1)=1.0;
cheb(1,2)=yC;
for jj=2:N-1;
cheb(1,jj+1)=2.0*yC*cheb(1,jj)-cheb(1,jj-1);
end
s=0;
for j=1:N
s = s + cheb(1,j)*real(V(j,58));% P:56, A:58, S:102
end
psi_hat(m)=s;
end

% plot
f1 = figure(1);
W = 4; H = 4;
set(f1,'PaperUnits','inches');set(f1,'PaperOrientation','portrait');
set(f1,'PaperSize',[H,W])    ;set(f1,'PaperPosition',[0,0,W,H]);
plot(real_c,imag_c,'ko',[0 1],[0 0]);axis([0 1 -1 0.1])
xlabel('c_{r}');ylabel('c_{i}');
print(f1,'-deps','-color','eig_vals.eps');

f2 = figure(2);
W = 4; H = 4;
set(f2,'PaperUnits','inches');set(f2,'PaperOrientation','portrait');
set(f2,'PaperSize',[H,W])    ;set(f2,'PaperPosition',[0,0,W,H]);
plot(y,abs(psi_hat),'k');% axis([0 1 -1 0.1])
xlabel('y');ylabel('psi');
print(f2,'-deps','-color','eig_modes.eps');

