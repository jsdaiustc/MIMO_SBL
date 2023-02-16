function [X,crb_DOD,crb_DOA]=signal(Mt,Mr,theta, alpha, SNR, K)
N_alpha=length(alpha);
for i1=1:N_alpha
    At=exp(-1i*pi*(0:Mt-1)'*sin(theta(i1)*pi/180));
    Ar=exp(-1i*pi*(0:Mr-1)'*sin(alpha(i1)*pi/180));
    A(:,i1)=kron(At,Ar);
    
%     Bt(:,i1)=kron(-1i*pi*(0:Mt-1)'*cos(theta(i1)*pi/180).*At,Ar);
%     Br(:,i1)=kron(At,-1i*pi*(0:Mr-1)'*cos(alpha(i1)*pi/180 ).*Ar);
    
Bt(:,i1)=kron((pi/180)*(-1i)*pi*cos(theta(i1)*pi/180)*((0:Mt-1)'.*At),Ar);
Br(:,i1)=kron(At,(pi/180)*(-1i)*pi*cos(alpha(i1)*pi/180)*((0:Mr-1)'.*Ar));
   
end
Vj=diag(sqrt((   (10*ones(N_alpha,1)).^(SNR/10)   )/2));
S=Vj*(randn(N_alpha,K)+1i*randn(N_alpha,K));
noise=sqrt(1/2)*(randn(Mt*Mr,K)+1i*randn(Mt*Mr,K));
X=A*S+noise;


% B=-1i*pi*(0:Mt-1)'*cos(theta*pi/180).*A;
 PP = (S*S'/K).';
 I =eye(Mt*Mr);
 F = I-A*inv( A'*A )*A';
 Ht = Bt'*F*Bt;
  Hr = Br'*F*Br;
 crb_DOD = mean( sqrt( diag(1/(2*K)*inv( real(Ht.*PP) ) ) ) );
 crb_DOA = mean( sqrt( diag(1/(2*K)*inv( real(Hr.*PP) ) ) ) );
