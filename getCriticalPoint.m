clear;
clc;
close all;

e=exp(1);
X(1,1)=rand()+2;
X(2,1)=rand()+2.5;
X(3,1)=rand()+4;
X(4,1)=rand()+5;
X(5,1)=rand()+6;
X(6,1)=rand()+6;
delta_t=1;
L=3100;        %????????
N=round(L/delta_t);  %????????
ts=zeros(N,1);
sigma=0.5;

p(1)=-0.04;
p(2)=-0.035;
p(3)=-0.03;
p(4)=-0.02;
p(5)=-0.01;
p(6)=0.006;
for i=7:27
    p(i)=(i-6)/50;
end
for i=1:27
    pr(i)=p(28-i);
end

T=50;
sample_num=10;
prob=zeros(25,1);
init_trans=[0.5,0.5;0.5,0.5];
for i=1:7
    init_emiss(1,i)=1/7;
    init_emiss(2,i)=1/7;
end
%%

D=[-1 0 0 1 0 0 0 0 0 0;...
   -1 -1 0 0 0 0 0 0 0 0 ;...
   1 0 1 0 0 0 0 0 0 0 ;...
   -1 0 0 -1 0 0 0 0 0 0;...
   0 -1 1 0 -1 1 0 0 0 0;...
   0 0 0 0 1 -1 1 0 0 0;...
   0 0 0 0 0 1 -1 0 0 0;...
   0 0 0 0 1 0 0 -1 0 1;...
   0 0 0 0 1 0 0 0 -1 0;...
   0 0 0 0 1 0 0 0 0 -1];

%%

for bb=1:T
    CC=zeros(6,27,sample_num);
    Changed=zeros(sample_num,26);
    for l=1:27
        q(l)=0.96^(1/abs(pr(l)));
        E=[-3/10*q(l) 0 0 0 0 0 0 0 0 0;...
            0 -5/10 0 0 0 0 0 0 0 0;...
            0 0 -3/5 0 0 0 0 0 0 0 ;...
            0 0 0 -4/5 0 0 0 0 0 0 ;...
            0 0 0 0 -5/5 0 0 0 0 0 ;...
            0 0 0 0 0 -6/5 0 0 0 0 ;...
            0 0 0 0 0 0 -7/5 0 0 0 ;...
            0 0 0 0 0 0 0 -8/5 0 0 ;...
            0 0 0 0 0 0 0 0 -9/5 0 ;...
            0 0 0 0 0 0 0 0 0 -10/5];
        J=D*E*inv(D);     
        J=[-2*q(l)/5 1-2*q(l)/5 1-2*q(l)/5 0 0 0; 2/5-q(l)/5 -q(l)/5-2/5 2/5-q(l)/5 0 0 0; q(l)/5-2/5 q(l)/5-3/5 q(l)/5-7/5 0 0 0; 0 -1/10 -1/10 -8/5 3/10 -3/10; 0 0 0 0 -21/10 1/10; 0 0 0 0 1/10 -21/10];
        for k=1:sample_num
            for i=1:N-1
                ts(i+1)=ts(i)+delta_t;
                eJ=e^(J*delta_t);
                for jj=1:10
                    X(jj,i+1)=eJ(jj,:)*X(:,i)+sigma*normrnd(0,1)*delta_t;
                end
            end
            CC(:,l,k)=X(:,3000);
        end
        if l>=2
            Changed(:,l-1)=changed_genes(CC,l);
        end
        if l>=3
            obs_seq=mat2cell(Changed(:,1:l-2),ones(1,sample_num))';
            [state_transi,emission]=hmmtrain(obs_seq,init_trans,init_emiss);
            aver_Pt=0;
            pi=[0.5,0.5];
            for i=1:sample_num
                aver_Pt=aver_Pt+pr_hmm(Changed(i,l-2:l-1),state_transi,...
                    emission,pi)/sample_num;
            end
            prob(l-2)=prob(l-2)+aver_Pt;
        end
                
    end    
    bb
end
prob=prob/bb;








