clc
clear all
%获取已知数据
m = 50000;%摄影比例尺分母
x0 = 0.0; %读入内方位元素
y0 = 0.0; 
fk = 153.24 / 1000.0; %单位为mm
DMZ = xlsread('地面坐标.xlsx','Sheet1','C4:E7');
XPZ = xlsread('像片坐标.xlsx','Sheet1','C4:D7');
XPZ = XPZ / 1000.0;   %单位转换,mm转为m 
n = size(DMZ,1); %控制点个数n 
%确定未知数的初始值
Xs = sum(DMZ(:,1)) / n; %Xs、Ys的初始值取控制点坐标的均值
Ys = sum(DMZ(:,2)) / n; 
Zs = sum(DMZ(:,3)) / n + m * fk; %Zs0=H=mf
t = 0.0; %角元素初始值为0
w = 0.0; 
k = 0.0; 
A = zeros(n * 2, 6); %A--误差方程式中的系数项(偏导数) 
L = zeros(n * 2, 1); %L--误差方程式中的常数项 

while (1) %迭代计算
%计算旋转矩阵R，ai、bi、ci为方向余旋
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
%逐点计算像片坐标近似值
%逐点计算误差方程式中的常数项和系数项并组成矩阵 
    for i = 1:n
%求像片近似坐标(x),(y)
    	xj = -fk * (a1 * (DMZ(i,1) - Xs) + b1 * (DMZ(i, 2) - Ys) + c1 * (DMZ(i,3) - Zs)) / (a3 * (DMZ(i, 1) - Xs) + b3 * (DMZ(i, 2) - Ys) + c3 * (DMZ(i, 3) - Zs));         
    	yj = -fk * (a2 * (DMZ(i,1) - Xs) + b2 * (DMZ(i, 2) - Ys) + c2 * (DMZ(i,3) - Zs)) / (a3 * (DMZ(i, 1) - Xs) + b3 * (DMZ(i, 2) - Ys) + c3 * (DMZ(i, 3) - Zs)); 
    	%A矩阵        
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
    	%L矩阵   
    	L(i * 2 - 1, 1) = XPZ(i, 1) - xj;          
    	L(i * 2, 1) = XPZ(i, 2) - yj;     
	  end 
	%解外方位元素改正值
	%X是一个列向量,存放一次计算之后外方位元素的改正数     
 	X = (A' * A) \ A' * L;  
	%修正外方位元素 
	Xs = X(1) + Xs;     
	Ys = X(2) + Ys;     
	Zs = X(3) + Zs;    
	t = X(4) + t;     
	w = X(5) + w;    
	k = X(6) + k; 
	%判断改正值是否小于给定限值
	%小于给定限值则退出循环 
	%若不满足要求则继续计算,直到满足要求为止    
	if abs(X(4)) < 0.001 && abs(X(5)) < 0.001 && abs(X(6)) < 0.001         
	 	break;     
	end 
end 
%输出所求得的参数
fprintf('Xs = %f\nYs = %f\nZs = %f\n', Xs, Ys, Zs); 
fprintf('t = %f\nw = %f\nk = %f\n', t, w, k);
fprintf('R = \n%f %f %f\n%f %f %f\n%f %f %f',a1,a2,a3,b1,b2,b3,c1,c2,c3);