function [DOD,DOA]=MIMO_SBL(Mt,Mr,Y,N_signal)

% Sparse Bayesian Approach for DOD and DOA Estimation With Bistatic MIMO Radar
% IEEE ACCESS 2020
% written by Zheng Cao

[M,L]=size(Y);
etc=10;


%% create the grid points according equation (8)
N_point=200; 
search_area_doa=(-90:180/N_point:90-180/N_point);
N_phi=floor(sqrt(N_point));
length_phi=length(search_area_doa);
N_un=1+floor(length_phi/N_phi-0.00001);
N_res=N_un*N_phi;
temp_add=(-90:180/N_res:90-180/N_res);
phi_add(1:N_res)=reshape( reshape(temp_add,N_un,N_phi)',N_res,1);
search_area_dod=phi_add(1:N_point);


%% Taylor expansion equation (4)
for i=1:N_point
    A_dod=exp(-1i*pi*(0:Mt-1)'*sin(search_area_dod(i)*pi/180));
    A_doa=exp(-1i*pi*(0:Mr-1)'*sin(search_area_doa(i)*pi/180));
    A(:,i)=kron(A_dod,A_doa);
    B_dod=-1i*pi*(0:Mt-1)'*cos((search_area_dod(i)*pi/180)).*A_dod;
    B_doa=-1i*pi*(0:Mr-1)'*cos((search_area_doa(i)*pi/180)).*A_doa;
    B1(:,i)=kron(B_dod,A_doa);
    B2(:,i)=kron(A_dod,B_doa);
end

%% initialization
a=0.0001;b=0.0001;
maxiter=200;
tol=1e-7;
alpha0=1;
delta_inv=mean(abs(A'*Y), 2)/(norm(A))^2;
converged = false;
iter = 0;


%% calculate mu and Sigma according to equation (20)
while ~converged
    iter = iter + 1;
    delta_last = delta_inv;  
    Phi=A;
    Phi_delta = Phi *  diag(delta_inv);
    V_temp= 1/alpha0*eye(M) +Phi_delta * Phi';  
    Sigma = diag(delta_inv) -Phi_delta' * (V_temp \Phi_delta);
    mu = alpha0 * Sigma * Phi' * Y;
    gamma1 = 1 - real(diag(Sigma)) ./ (delta_inv);
    
    
    %% update signal precision according to equation (22)
    temp=sum( mu.*conj(mu),2) + L*real(diag(Sigma));
    delta_inv= ( b+ real(temp)  )/(L+a);
    
    
    
    %% update noise precision according to equation (21)
    resid=Y-Phi*mu;
    alpha0=( L*M+ (1))/( b + norm(resid, 'fro')^2  +  (L / alpha0) * sum(gamma1)   );
    
    
    
    %% stopping criteria
    erro=norm(delta_inv - delta_last)/norm(delta_last);
    if ((erro < tol)&&iter>100) || (iter >= maxiter)
        converged = true;
    end

    
    %% off-grid refining
    Pm=mean(mu.*conj(mu),2);
    [~,sort_ind]=sort(Pm, 'descend');
    idx=sort_ind(1:etc);
    etc_length=length(idx);    
    dds=180/N_phi;
    if iter<=50
        resol1=dds;
    else
        resol1=1;    
    end
    
    % fix DOD grid gap according to (23)
    [temp1]=update_beta(idx,A,L,Y,etc_length,B1,mu,Sigma);
    search_area_dod(idx)=search_area_dod(idx)+sign(temp1')*resol1/100+temp1'*180/pi;
      
    % fix DOA grid gap according to (23)
    [temp2]=update_beta(idx,A,L,Y,etc_length,B2,mu,Sigma);
    search_area_doa(idx)=search_area_doa(idx)+   sign(temp2')*resol1/100 + temp2'*180/pi;
 
    temp_dod=search_area_dod(idx);
    temp_doa=search_area_doa(idx);
        
    for i=1:etc
        A_dod=exp(-1i*pi*(0:Mt-1)'*sin(temp_dod(i)*pi/180));
        A_doa=exp(-1i*pi*(0:Mr-1)'*sin(temp_doa(i)*pi/180));
        A(:,sort_ind(i))=kron(A_dod,A_doa);
        B_dod=-1i*pi*(0:Mt-1)'*cos((temp_dod(i)*pi/180)).*A_dod;
        B_doa=-1i*pi*(0:Mr-1)'*cos((temp_doa(i)*pi/180)).*A_doa;
        B1(:,sort_ind(i))=kron(B_dod,A_doa);
        B2(:,sort_ind(i))=kron(A_dod,B_doa);
    end
    
    
    %% remove the fause peaks beside the true DOD/DOA peaks 
    if iter>=101 && iter<102       
        ind_temp=1;
        kk=1;
        while length(ind_temp)<N_signal
            if min(abs(temp_doa(ind_temp)-temp_doa(kk+1)))>3 &  min(abs(temp_dod(ind_temp)-temp_dod(kk+1)))>3 
                ind_temp=[ind_temp,kk+1];
                kk=kk+1;
            else
                kk=kk+1;
            end
        end        
        idx_2=idx(ind_temp);
        etc=N_signal;
        search_area_dod=search_area_dod(idx_2);
        search_area_doa=search_area_doa(idx_2);
        A=A(:,idx_2);
        B1=B1(:,idx_2);
        B2=B2(:,idx_2);
        delta_inv=delta_inv(idx_2)*0 + 1;
        alpha0=1;
    end      
end

DOD1=search_area_dod;
DOA1=search_area_doa;
[DOD,J]=sort(DOD1);
DOA=DOA1(J);





