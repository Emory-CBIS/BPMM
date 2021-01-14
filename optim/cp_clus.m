function cp = cp_clus(x,percen,maxdist)
T = size(x,2);
    cp1 = x;
    cp_cut = round(sum(x)/50);
    cp2 = quantile(cp1,percen);
    cp1(cp1<=cp2) = 0;
    tbl_jum = zeros(T,2);
    tbl_jum(:,2) = cp1;
    tbl_jum(:,1) = 1:T;
tbl_jum2=tbl_jum(tbl_jum(:,2)~=0,:); 
dist_max = maxdist;

% tbl_jum2=tbl_jum2(tbl_jum2(:,2)~=1,:);
n=T;
% jum_pop2/3 grouping results
jum_pop=zeros(size(tbl_jum2,1),size(tbl_jum2,1)); %the grouping of all possible change points
jum_pop3=zeros(size(tbl_jum2,1),2); %the value of each group and sum of all prob
jum_prob=zeros(1,size(tbl_jum2,1));
if tbl_jum2(1,1) ~= n
for i=1:size(tbl_jum2,1)
   if size(tbl_jum2,1)==1
            jum_pop(1,1)=tbl_jum2(1,1);
            jum_pop3(1,1)=tbl_jum2(i,1);
            jum_pop3(1,2)=tbl_jum2(i,2);
   else 
   if i==1 
       n_row=1;
       n_col=1;
       jum_pop(1,1)=tbl_jum2(1,1);
       jum_prob(1,1)=tbl_jum2(1,2);
   elseif i<size(tbl_jum2,1)
       if tbl_jum2(i,1)-tbl_jum2(i-1,1)<=dist_max
           n_col=n_col+1;
           jum_pop(n_row,n_col)=tbl_jum2(i,1);
           jum_prob(1,n_col)=tbl_jum2(i,2);
       else
           %[~,ind]=max(jum_prob);
           jum_prob_sum = sum(jum_prob);
           jum_ind = round(jum_prob_sum/2); %the value of the group is the median one
           ind = 0;
           inds = 0;
           while inds<jum_ind
               ind = ind +1;
               inds = inds + jum_prob(ind);
           end
           jum_pop3(n_row,1)=jum_pop(n_row,ind);
           %[ind,~]=max(jum_prob);
           %jum_pop3(n_row,1)=round(mean(jum_pop(n_row,jum_prob(1,:)==ind)));
           %jum_pop3(n_row,1)=round(jum_pop(n_row,:)*jum_prob'/sum(jum_prob));
           jum_pop3(n_row,2)=sum(jum_prob);
           jum_prob=zeros(1,size(tbl_jum2,1));
           n_col=1;
           n_row=n_row+1;
           jum_pop(n_row,n_col)=tbl_jum2(i,1);
           jum_prob(1,n_col)=tbl_jum2(i,2);
       end
   else
       if tbl_jum2(i,1)-tbl_jum2(i-1,1)<=dist_max
           n_col=n_col+1;
           jum_pop(n_row,n_col)=tbl_jum2(i,1);
           jum_prob(1,n_col)=tbl_jum2(i,2);
           %[~,ind]=max(jum_prob);
           jum_prob_sum = sum(jum_prob);
           jum_ind = round(jum_prob_sum/2);
           ind = 0;
           inds = 0;
           while inds<jum_ind
               ind = ind +1;
               inds = inds + jum_prob(ind);
           end
           jum_pop3(n_row,1)=jum_pop(n_row,ind);
           jum_pop3(n_row,2)=sum(jum_prob);
       else
           %[~,ind]=max(jum_prob);
           jum_prob_sum = sum(jum_prob);
           jum_ind = round(jum_prob_sum/2);
           ind = 0;
           inds = 0;
           while inds<jum_ind
               ind = ind +1;
               inds = inds + jum_prob(ind);
           end
           jum_pop3(n_row,1)=jum_pop(n_row,ind);
           n_col=1;
           n_row=n_row+1;
           jum_pop3(n_row,1)=tbl_jum2(i,1);
           jum_pop3(n_row,2)=tbl_jum2(i,2);
           jum_pop(n_row,1)=tbl_jum2(i,1);
       end
   end
    end
end
end

jum_pop2 = jum_pop(all(jum_pop==0,2)==0,:);
n_col=size(tbl_jum2,1)-min(sum(jum_pop2==0,2));
jum_pop2=jum_pop2(:,1:n_col);
for i=1:size(jum_pop2,1)
jum_pop3(i,2)=sum(tbl_jum2(ismember(tbl_jum2(:,1),jum_pop2(i,:)),2));
end
jum_pop3= jum_pop3(all(jum_pop3==0,2)==0,:);
[~,I]=sort(jum_pop3(:,2),'descend');
jum_pop4 = jum_pop3(I,:);
cp = sort(jum_pop4(jum_pop4(:,2)>cp_cut,1)');
 
    
end


