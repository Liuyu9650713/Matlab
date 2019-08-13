function [data,flag,active] = ID3_Handle(  )
%% ID3�㷨����Ԥ�������ַ���ת��Ϊ0,1����

% ���������
% data�� ת�����0,1����
% flag: ���Ժ�Label��
% active: ����������ȫ1��

%% ��ȡ����
text={'���' '����' '�й���' '���Լ��ķ���' '�Ŵ����' '���'
        '1'  '����'   '��'       '��'        'һ��'   '��'
        '2'  '����'   '��'       '��'         '��'    '��'
        '3'  '����'   '��'       '��'         '��'    '��'
        '4'  '����'   '��'       '��'        'һ��'   '��'
        '5'  '����'   '��'       '��'        'һ��'   '��'
        '6'  '����'   '��'       '��'        'һ��'   '��'
        '7'  '����'   '��'       '��'         '��'    '��'
        '8'  '����'   '��'       '��'         '��'    '��'
        '9'  '����'   '��'       '��'       '�ǳ���'   '��'
        '10' '����'   '��'       '��'       '�ǳ���'   '��'
        '11' '����'   '��'       '��'       '�ǳ���'   '��'
        '12' '����'   '��'       '��'         '��'     '��'
        '13' '����'   '��'       '��'         '��'     '��'
        '14' '����'   '��'       '��'       '�ǳ���'   '��'
        '15' '����'   '��'       '��'        'һ��'    '��'}
    
    matrix=text(2:end,2:end);
    flag=text(1,2:end-1);
    active=ones(1,length(flag));
    
%% ����ת��
 [rows,cols]=size(matrix);
 data=zeros(rows,cols);
 for i=1:rows
     data(i,:)=cellfun(@translate,matrix(i,:));
 end
end
 

function tag= translate(matrix)
    if strcmp(matrix,'����') ||strcmp(matrix,'��')...
        ||strcmp(matrix,'һ��')
        tag =0;
        return ;
    end
     if strcmp(matrix,'����') ||strcmp(matrix,'��')...
        ||strcmp(matrix,'��')
        tag =1;
        return ;
    end
    tag =2;
end
 
