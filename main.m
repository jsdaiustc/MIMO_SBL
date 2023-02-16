clear;

Mt=6;
Mr=6;
N_snap=20;
SNR=20;
DOD_real=[-17.4, 12.7 ]; 
DOA_real=[-6.5,   20.2]; 

[DOD_real,J]=sort(DOD_real);
DOA_real=DOA_real(J);

N_signal=length(DOD_real);
Y=signal(Mt,Mr,DOD_real,DOA_real,SNR, N_snap);
               
[DOD,DOA]=MIMO_SBL(Mt,Mr,Y,N_signal)
