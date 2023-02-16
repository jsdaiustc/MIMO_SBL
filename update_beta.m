function [temp1,idx]=update_beta(idx,A,T,Y,k,B,mu,Sigma)

    BHB = B(:,idx)' * B(:,idx);  
    
    P1 = real( conj(BHB) .* ((mu(idx,:) * mu(idx,:)') + T * Sigma(idx,idx)) );%conjÇó¸´¹²éî
%     A1 = A ;%+ B2*diag(beta2);
    
    
        v1 = sum( real(conj(mu(idx,:)) .* (B(:,idx)' * (Y - A * mu))),2) ;   
    
    v1 = v1 - T * real(diag(B(:,idx)' * A * Sigma(:,idx)));
    
 
    
%     eigP=svd(P1);
%     if eigP(end)/eigP(1)>1e-5
           temp1 =  P1 \ v1;
%      else
%            temp1=v1./diag(P1);
%     end
     
%    if iter<200 
%       theld=0.01       
%    else
%        theld=0.01
%    end

%      theld=0.005;   
%      detheta=temp1'*180/pi;
%      ind_sma=find(abs(detheta)<theld);
%      detheta(ind_sma)=detheta(ind_sma)./abs(detheta(ind_sma))*0.005;
%      temp1= detheta'*pi/180;

    