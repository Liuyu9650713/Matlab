function [data,flag,active] = ID3_Handle(  )
%% ID3算法数据预处理，把字符串转换为0,1编码

% 输出参数：
% data： 转换后的0,1矩阵；
% flag: 属性和Label；
% active: 属性向量，全1；

%% 读取数据
text={'序号' '年龄' '有工作' '有自己的房子' '信贷情况' '类别'
        '1'  '青年'   '否'       '否'        '一般'   '否'
        '2'  '青年'   '否'       '否'         '好'    '否'
        '3'  '青年'   '是'       '否'         '好'    '是'
        '4'  '青年'   '是'       '是'        '一般'   '是'
        '5'  '青年'   '否'       '否'        '一般'   '否'
        '6'  '中年'   '否'       '否'        '一般'   '否'
        '7'  '中年'   '否'       '否'         '好'    '否'
        '8'  '中年'   '是'       '是'         '好'    '是'
        '9'  '中年'   '否'       '是'       '非常好'   '是'
        '10' '中年'   '否'       '是'       '非常好'   '是'
        '11' '老年'   '否'       '是'       '非常好'   '是'
        '12' '老年'   '否'       '是'         '好'     '是'
        '13' '老年'   '是'       '否'         '好'     '是'
        '14' '老年'   '是'       '否'       '非常好'   '是'
        '15' '老年'   '否'       '否'        '一般'    '否'}
    
    matrix=text(2:end,2:end);
    flag=text(1,2:end-1);
    active=ones(1,length(flag));
    
%% 数据转换
 [rows,cols]=size(matrix);
 data=zeros(rows,cols);
 for i=1:rows
     data(i,:)=cellfun(@translate,matrix(i,:));
 end
end
 

function tag= translate(matrix)
    if strcmp(matrix,'青年') ||strcmp(matrix,'否')...
        ||strcmp(matrix,'一般')
        tag =0;
        return ;
    end
     if strcmp(matrix,'中年') ||strcmp(matrix,'好')...
        ||strcmp(matrix,'是')
        tag =1;
        return ;
    end
    tag =2;
end
 
