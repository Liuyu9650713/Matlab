function ID3_Printtree(tree)
global nodeid;
nodeid=0;      %主要用来统计节点个数
if isempty(tree) 
    disp('空树！');
    return ;
end
queue = queue_push([],tree);   %将根节点输入到队列中

while ~isempty(queue) % 队列不为空
     [node,queue] = queue_pop(queue); % 出队列

     visit(node,length(queue));  %访问该结点  既将此节点的属性输出
     if ~strcmp(node.lchild,'NULL') % 左子树不为空
        queue = queue_push(queue,node.lchild); % 进队
     end
     if ~strcmp(node.rchild,'NULL') % 右子树不为空
        queue = queue_push(queue,node.rchild); % 进队
     end
     if ~strcmp(node.mchild,'NULL') % 中子树不为空
        queue = queue_push(queue,node.mchild); % 进队
     end
end
end


function visit(node,length)
%% 访问node 节点，并把其设置值为nodeid的节点
    global nodeid;
    if isleaf(node)
        nodeid=nodeid+1;
        fprintf('叶子节点，node: %d\t，属性值: %s\n', nodeid, node.name);
    else 
        nodeid=nodeid+1;
         if strcmp(node.mchild,'NULL')   %如果该节点属性有两种情况时为第一种输出方式，三种时为第二种输出方式
            fprintf('node: %d\t属性值: %s\t，左子树为节点：node%d，右子树为节点：node%d\n',nodeid, node.name,nodeid+length+1,nodeid+length+2);
        else
            fprintf('node: %d\t属性值: %s\t，左子树为节点：node%d，右子树为节点：node%d,中子树为节点：node%d\n',nodeid, node.value,nodeid+length+1,nodeid+length+2,nodeid+3);
         end
    end
end

function tag = isleaf(node)
%% 是否是叶子节点
    if strcmp(node.lchild,'NULL') && strcmp(node.rchild,'NULL') && strcmp(node.mchild,'NULL') % 左右中都为空
        tag =1;
    else
        tag=0;
    end
end

function  [ item,newqueue ]=queue_pop(queue)
%% 出队
if isempty(queue)
    disp('队列为空，不能访问！');
    return;
end
item = queue(1); % 第一个元素弹出
newqueue=queue(2:end); % 往后移动一个元素位置
end

function [ newqueue ] = queue_push( queue,item )
%% 进队
newqueue=[queue,item];

end