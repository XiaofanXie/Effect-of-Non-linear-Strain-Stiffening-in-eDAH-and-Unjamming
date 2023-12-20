%%springpot (fractional element) fit
%%Frank Sauer 17.03.2021

function[mu,alpha]= springpot(f,G1,G2);
%function
 model = @(x,xdata)(x(1)^(1-x(2)))*(i*2*pi*xdata).^(x(2));
 fun = @(x,xdata) [real(model(x,xdata)); imag(model(x,xdata))];
%starting value  
 x0=[10,0.5];
 lb=[0,0];
 ub=[10000,1];
 %
 xdata = f;
 ydata= [G1;G2];
%fit
 x = lsqcurvefit(fun,x0,xdata,ydata,lb,ub);
 mu=x(1);
 alpha=x(2);

% % test data
%  f = [5,10,20,50,100,200];
%  G1 = [42.2495105655658,49.3626884577083,65.2854872233429,92.8565087858124,133.199833550814,277.123784159478];
%  G2 = [37.6636433486276,61.4187347153510,101.310689370004,197.694190305441,379.603759751271,685.949615597812];

% %   plot test data
close all
times = linspace(xdata(1),xdata(end));
plot(xdata,(ydata),'o')
hold on;
plot(times,fun(x,times))
pause;
