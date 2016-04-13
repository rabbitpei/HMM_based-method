function [P]=pro(features,t,ii,jj)
%n features
%calculate pro(jj|jj), ie. P(s(t)=jj|s(t-1)=ii)
%P(s(t)=jj|s(t-1)=ii)=P(s(t-1)=ii,s(t)=jj)/P(s(t-1)=ii), ie. P(s1|s2)=P(s2,s1)/P(s2)
if nargin==3
    jj=0;
end

[n,tp,sample_num]=size(features);
for i=1:n
    vc=features(i,t-1,:);
    [mu(i),sigma(i)]=normfit(vc(:));
end
samples2=reshape(features(:,t-1,:),n,sample_num);
samples1=reshape(features(:,t,:),n,sample_num);
observed_pre=zeros(n,sample_num);
observed_post=zeros(n,sample_num);
for i=1:n
    for j=1:sample_num
        if samples2(i,j)>mu(i)-sigma(i)&&samples2(i,j)<mu(i)+sigma(i)
            observed_post(i,j)=1;
        end
        if samples1(i,j)>mu(i)-sigma(i)&&samples1(i,j)<mu(i)+sigma(i)
            observed_pre(i,j)=1;
        end
    end
end
coin1=0;
coin2=0;
for s=1:sample_num
    if sum(observed_post(:,s))==ii
        coin2=coin2+1;
        if sum(observed_pre(:,s))==jj
            coin1=coin1+1;
        end
    end
end

if nargin==3
    P=coin2/sample_num;
else 
    P=coin1/sample_num;
end
    





