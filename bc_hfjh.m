clc
clear all
%��ȡ��֪����
m = 50000;%��Ӱ�����߷�ĸ
x0 = 0.0; %�����ڷ�λԪ��
y0 = 0.0; 
fk = 153.24 / 1000.0; %��λΪmm
DMZ = xlsread('��������.xlsx','Sheet1','C4:E7');
XPZ = xlsread('��Ƭ����.xlsx','Sheet1','C4:D7');
XPZ = XPZ / 1000.0;   %��λת��,mmתΪm 
n = size(DMZ,1); %���Ƶ����n 
%ȷ��δ֪���ĳ�ʼֵ
Xs = sum(DMZ(:,1)) / n; %Xs��Ys�ĳ�ʼֵȡ���Ƶ�����ľ�ֵ
Ys = sum(DMZ(:,2)) / n; 
Zs = sum(DMZ(:,3)) / n + m * fk; %Zs0=H=mf
t = 0.0; %��Ԫ�س�ʼֵΪ0
w = 0.0; 
k = 0.0; 
A = zeros(n * 2, 6); %A--����ʽ�е�ϵ����(ƫ����) 
L = zeros(n * 2, 1); %L--����ʽ�еĳ����� 

while (1) %��������
%������ת����R��ai��bi��ciΪ��������
    a1 = cos(t) * cos(k) - sin(t) * sin(w) * sin(k);     
    a2 = -cos(t) * sin(k) - sin(t) * sin(w) * cos(k);    
    a3 = -sin(t) * cos(w);     
    b1 = cos(w) * sin(k);     
    b2 = cos(w) * cos(k);     
    b3 = -sin(w);     
    c1 = sin(t) * cos(k) + cos(t) * sin(w) * sin(k);     
    c2 = -sin(t) * sin(k) + cos(t) * sin(w) * cos(k);     
    c3 = cos(t) * cos(w);
    R= [a1 a2 a3;
        b1 b2 b3;
        c1 c2 c3];
%��������Ƭ�������ֵ
%����������ʽ�еĳ������ϵ�����ɾ��� 
    for i = 1:n
%����Ƭ��������(x),(y)
    	xj = -fk * (a1 * (DMZ(i,1) - Xs) + b1 * (DMZ(i, 2) - Ys) + c1 * (DMZ(i,3) - Zs)) / (a3 * (DMZ(i, 1) - Xs) + b3 * (DMZ(i, 2) - Ys) + c3 * (DMZ(i, 3) - Zs));         
    	yj = -fk * (a2 * (DMZ(i,1) - Xs) + b2 * (DMZ(i, 2) - Ys) + c2 * (DMZ(i,3) - Zs)) / (a3 * (DMZ(i, 1) - Xs) + b3 * (DMZ(i, 2) - Ys) + c3 * (DMZ(i, 3) - Zs)); 
    	%A����        
    	Zbar = a3 * (DMZ(i, 1) - Xs) + b3 * (DMZ(i, 2) - Ys) + c3 * (DMZ(i, 3) - Zs);        
    	A(i * 2 - 1, 1) = (a1 * fk + a3 * XPZ(i, 1)) / Zbar;        
    	A(i * 2 - 1, 2) = (b1 * fk + b3 * XPZ(i, 1)) / Zbar;         
    	A(i * 2 - 1, 3) = (c1 * fk + c3 * XPZ(i, 1)) / Zbar;        
    	A(i * 2 - 1, 4) = XPZ(i, 2) * sin(w) - (XPZ(i, 1) * (XPZ(i, 1) * cos(k) - XPZ(i, 2) * sin(k)) / fk + fk * cos(k)) * cos(w);        
    	A(i * 2 - 1, 5) = -fk * sin(k) - XPZ(i, 1) * (XPZ(i, 1) * sin(k) + XPZ(i, 2) * cos(k)) / fk;         
    	A(i * 2 - 1, 6) = XPZ(i, 2);        
    	A(i * 2, 1) = (a2 * fk + a3 * XPZ(i, 2)) / Zbar;         
    	A(i * 2, 2) = (b2 * fk + b3 * XPZ(i, 2)) / Zbar;         
    	A(i * 2, 3) = (c2 * fk + c3 * XPZ(i, 2)) / Zbar;        
    	A(i * 2, 4) = XPZ(i, 1) * sin(w) - (XPZ(i, 2) * (XPZ(i, 1) * cos(k) - XPZ(i, 2) * sin(k)) / fk - fk * sin(k)) * cos(w);         
    	A(i * 2, 5) = -fk * cos(k) - XPZ(i, 2) * (XPZ(i, 1) * sin(k) + XPZ(i, 2) * cos(k)) / fk;        
    	A(i * 2, 6) = -XPZ(i, 1);         
    	%L����   
    	L(i * 2 - 1, 1) = XPZ(i, 1) - xj;          
    	L(i * 2, 1) = XPZ(i, 2) - yj;     
	  end 
	%���ⷽλԪ�ظ���ֵ
	%X��һ��������,���һ�μ���֮���ⷽλԪ�صĸ�����     
 	X = (A' * A) \ A' * L;  
	%�����ⷽλԪ�� 
	Xs = X(1) + Xs;     
	Ys = X(2) + Ys;     
	Zs = X(3) + Zs;    
	t = X(4) + t;     
	w = X(5) + w;    
	k = X(6) + k; 
	%�жϸ���ֵ�Ƿ�С�ڸ�����ֵ
	%С�ڸ�����ֵ���˳�ѭ�� 
	%��������Ҫ�����������,ֱ������Ҫ��Ϊֹ    
	if abs(X(4)) < 0.001 && abs(X(5)) < 0.001 && abs(X(6)) < 0.001         
	 	break;     
	end 
end 
%�������õĲ���
fprintf('Xs = %f\nYs = %f\nZs = %f\n', Xs, Ys, Zs); 
fprintf('t = %f\nw = %f\nk = %f\n', t, w, k);
fprintf('R = \n%f %f %f\n%f %f %f\n%f %f %f',a1,a2,a3,b1,b2,b3,c1,c2,c3);