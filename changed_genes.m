function changed=changed_genes(features,t)
[n,tp,sample_num]=size(features);
for i=1:n
    vc=features(i,t-1,:);
    [mu(i),sigma(i)]=normfit(vc(:));
end
samples=reshape(features(:,t,:),n,sample_num);
observed=zeros(1,sample_num);
for j=1:sample_num
    for i=1:n
        if samples(i,j)>mu(i)-sigma(i)&&samples(i,j)<mu(i)+sigma(i)
            observed(j)=observed(j)+1;
        end
    end
end
changed=observed+1;