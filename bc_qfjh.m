clc
clear all
%���������������
Zxy=xlsread('��Ƭ.xlsx','Sheet1','C4:D7');
Yxy=xlsread('��Ƭ.xlsx','Sheet1','C4:D7');
n = size(Zxy,1); %�����������n 
XXZP=[n,3];
%�����ڷ�λԪ�ء��ⷽλԪ��
f=150.0/1000.0;
x0=0.0;
y0=0.0;
%��Ƭ�ⷽλԪ��
t1=0.000215;  w1=0.029064;  k1=0.095247;
Xs1=49997.70;  Ys1=49997.29;  Zs1=20000.02;
%��Ƭ�ⷽλԪ��
t2=0.014434;  w2=0.046018;  k2=0.110469;                                   
Xs2=58968.29;  Ys2=50702.44;  Zs2=20304.43;

%ȷ����������������R����
%��Ƭ
a11 = cos(t1) * cos(k1) - sin(t1) * sin(w1) * sin(k1);     
a21 = -cos(t1) * sin(k1) - sin(t1) * sin(w1) * cos(k1);    
a31 = -sin(t1) * cos(w1);     
b11 = cos(w1) * sin(k1);
b21 = cos(w1) * cos(k1);     
b31 = -sin(w1);
c11 = sin(t1) * cos(k1) + cos(t1) * sin(w1) * sin(k1);     
c21 = -sin(t1) * sin(k1) + cos(t1) * sin(w1) * cos(k1);     
c31 = cos(t1) * cos(w1);
R1= [a11 a21 a31;
     b11 b21 b31;
     c11 c21 c31];

%��Ƭ
a12 = cos(t2) * cos(k2) - sin(t2) * sin(w2) * sin(k2);     
a22 = -cos(t2) * sin(k2) - sin(t2) * sin(w2) * cos(k2);    
a32 = -sin(t2) * cos(w2);     
b12 = cos(w2) * sin(k2);     
b22 = cos(w2) * cos(k2);     
b32 = -sin(w2);     
c12 = sin(t2) * cos(k2) + cos(t2) * sin(w2) * sin(k2);     
c22 = -sin(t2) * sin(k2) + cos(t2) * sin(w2) * cos(k2);     
c32 = cos(t2) * cos(w2);
R2= [a12 a22 a32;
     b12 b22 b32;
     c12 c22 c32];
 
%������߷���
Bx = Xs2 - Xs1;
By = Ys2 - Ys1;
Bz = Zs2 - Zs1;


for i=1:n %�ֱ������������
    %������ռ丨������
    l1=[Zxy(i,1);Zxy(i,2);-f];
    l2=[Yxy(i,1);Yxy(i,2);-f];
    FZ=R1*l1;
    FY=R2*l2;
    %����ͶӰϵ��
    N1 = (Bx*FY(3,1) - Bz*FY(1,1)) / (FZ(1,1)*FY(3,1) - FZ(3,1)*FY(1,1));
	N2 = (Bx*FZ(3,1) - Bz*FZ(1,1)) / (FZ(1,1)*FY(3,1) - FZ(3,1)*FY(1,1));
    %�������������ϵ����
    dX = N1*FZ(1,1);
	dY = 0.5*(N1*FZ(2,1) + N2*FY(2,1) + By);
	dZ = N1*FZ(3,1);
    %��������ĵ�������
    XXZP(i,1)= Xs1 + dX;
	XXZP(i,2) = Ys1 + dY;
	XXZP(i,3) = Zs1 + dZ;
end

fprintf('�������� = \n');
fprintf('%f %f %f\n', XXZP(1,1), XXZP(1,2), XXZP(1,3));
fprintf('%f %f %f\n', XXZP(2,1), XXZP(2,2), XXZP(2,3));
fprintf('%f %f %f\n', XXZP(3,1), XXZP(3,2), XXZP(3,3));
fprintf('%f %f %f\n', XXZP(4,1), XXZP(4,2), XXZP(4,3));
