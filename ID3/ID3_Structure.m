function [tree]=ID3_Structure(data,flag,active)

%% 检查数据及定义树
if(isempty(data))
    error('必须提供数据！');
end

%树的结构定义
tree=struct('name','NULL','rchild','NULL','lchild','NULL','mchild','NULL');

volume=length(data(:,1));  %样本容量
numberactive= length(active);  

lastcolum=sum(data(:,end)); %类别中‘是’出现的总量

%% 特殊节点情况

%第一种 类别全部为‘是’的情况
if lastcolum==volume
    tree.name='true';
    return
end

%第一种 类别全部为‘否’的情况
if lastcolum==0
    tree.name='false';
    return
end

%第三种 无属性可继续划分的情况
if sum(active)==0
    if(lastcolum>=volume/2)
        tree.name='true';
    else
        tree.name='false';
    end
    return
end

%% 经验熵
p1=lastcolum/volume;
if p1==0
    H_p1=0;
else
    H_p1=-1*p1*log2(p1);
end
p2=(volume-lastcolum)/volume;
if p2==0
    H_p2=0;
else
    H_p2=-1*p2*log2(p2);
end
H_D=H_p1+H_p2;

%% 计算信息增益并选取最大信息增益的属性

gains=-1*ones(1,numberactive);

for j=1:numberactive
    if active(j)==0
        continue
    end
     a=tabulate(data(:,j));
     b=a(:,1);
     c=length(b);
     s=zeros(1,c);
     s_C=zeros(1,c);
     for i=1:volume
         for k=1:c
             if data(i,j)==b(k)
                 s(k)=s(k)+1;
                 if data(i,end)
                     s_C(k)=s_C(k)+1;
                 end
             end
         end
     end
     
     H=zeros(1,c);
     for k=1:c
         if s(k)==0
             p1=0;
         else
             p1=s_C(k)/s(k);
         end
         if p1==0
             eq_q1=0;
         else
             eq_q1=-1*(p1)*log2(p1);
         end
         if s(k)==0
             p2=0;
         else
             p2=(s(k)-s_C(k))/s(k);
         end
         if p2==0;
             eq_q2=0;
         else
             eq_q2=-1*(p2)*log2(p2);
         end
         H(k)=(eq_q1+eq_q2)*s(k)/volume;
     end
     gains(j)=H_D-sum(H);                            
end
% 选出最大熵的位置；
[~, local] = max(gains);

% 设置相应值
tree.name= flag{local};
% 去活跃状态
active(local) = 0;


% 根据local把数据进行分组
a=tabulate(data(:,local));
b=a(:,1);
c=length(b);
for k=1:c
    if k==1
        if (isempty(data(data(:,local)==b(k),:)));
               leaf = struct('name', 'NULL', 'lchild', 'NULL', 'right', 'NULL','mchild','NULL');
           if (lastcolum>=volume/ 2);
               leaf.name = 'true';
           else
               leaf.name = 'false';
           end
          tree.lchild = leaf;
       else
         % 递归
             tree.lchild = ID3_Structure(data(data(:,local)==b(k),:), flag,active);
        end 
    end

    if k==2
        if (isempty(data(data(:,local)==b(k),:)));
               leaf = struct('name', 'NULL', 'lchild', 'NULL', 'rchild', 'NULL','mchild','NULL');
           if (lastcolum>=volume/ 2);
               leaf.name = 'true';
           else
               leaf.name = 'false';
           end
          tree.lchild = leaf;
       else
         % 递归
             tree.rchild = ID3_Structure(data(data(:,local)==b(k),:), flag,active);
        end 
    end
    if k==3
        if (isempty(data(data(:,local)==b(k),:)));
               leaf = struct('name', 'NULL', 'lchild', 'NULL', 'rchild', 'NULL','mchild','NULL');
           if (lastcolum>=volume/ 2);
               leaf.name = 'true';
           else
               leaf.name = 'false';
           end
          tree.lchild = leaf;
       else
         % 递归
             tree.mchild = ID3_Structure(data(data(:,local)==b(k),:), flag,active);
        end 
    end
end
return
end