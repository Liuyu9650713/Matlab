function ID3_Printtree(tree)
global nodeid;
nodeid=0;      %��Ҫ����ͳ�ƽڵ����
if isempty(tree) 
    disp('������');
    return ;
end
queue = queue_push([],tree);   %�����ڵ����뵽������

while ~isempty(queue) % ���в�Ϊ��
     [node,queue] = queue_pop(queue); % ������

     visit(node,length(queue));  %���ʸý��  �Ƚ��˽ڵ���������
     if ~strcmp(node.lchild,'NULL') % ��������Ϊ��
        queue = queue_push(queue,node.lchild); % ����
     end
     if ~strcmp(node.rchild,'NULL') % ��������Ϊ��
        queue = queue_push(queue,node.rchild); % ����
     end
     if ~strcmp(node.mchild,'NULL') % ��������Ϊ��
        queue = queue_push(queue,node.mchild); % ����
     end
end
end


function visit(node,length)
%% ����node �ڵ㣬����������ֵΪnodeid�Ľڵ�
    global nodeid;
    if isleaf(node)
        nodeid=nodeid+1;
        fprintf('Ҷ�ӽڵ㣬node: %d\t������ֵ: %s\n', nodeid, node.name);
    else 
        nodeid=nodeid+1;
         if strcmp(node.mchild,'NULL')   %����ýڵ��������������ʱΪ��һ�������ʽ������ʱΪ�ڶ��������ʽ
            fprintf('node: %d\t����ֵ: %s\t��������Ϊ�ڵ㣺node%d��������Ϊ�ڵ㣺node%d\n',nodeid, node.name,nodeid+length+1,nodeid+length+2);
        else
            fprintf('node: %d\t����ֵ: %s\t��������Ϊ�ڵ㣺node%d��������Ϊ�ڵ㣺node%d,������Ϊ�ڵ㣺node%d\n',nodeid, node.value,nodeid+length+1,nodeid+length+2,nodeid+3);
         end
    end
end

function tag = isleaf(node)
%% �Ƿ���Ҷ�ӽڵ�
    if strcmp(node.lchild,'NULL') && strcmp(node.rchild,'NULL') && strcmp(node.mchild,'NULL') % �����ж�Ϊ��
        tag =1;
    else
        tag=0;
    end
end

function  [ item,newqueue ]=queue_pop(queue)
%% ����
if isempty(queue)
    disp('����Ϊ�գ����ܷ��ʣ�');
    return;
end
item = queue(1); % ��һ��Ԫ�ص���
newqueue=queue(2:end); % �����ƶ�һ��Ԫ��λ��
end

function [ newqueue ] = queue_push( queue,item )
%% ����
newqueue=[queue,item];

end