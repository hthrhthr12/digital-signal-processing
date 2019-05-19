%%% Initialization
clear variables; clc
close all
c=.7*2*(1+.95^2-2*.95*cos(pi/3))/(17-8*cos(pi/8));%normalization
p1=.95*exp(1j*pi/3); p2=0.3;%poles
z1=4*exp(1j*pi/8); z2=0.5;%zeros
%%%%%%% Define 8 stable, real and casual systems
%%% ROC={z: abs(z)>0.95}
H{1}=@(z)c*((1-z2*z.^-1).*(1-z1*z.^-1).*(1-z1'*z.^-1))./...
    ((1-p2*z.^-1).*(1-p1*z.^-1).*(1-p1'*z.^-1));
H{2}=@(z)-H{1}(z);
% zplane(z,p)
H{3}=@(z)c*((z.^-1-z2).*(1-z1*z.^-1).*(1-z1'*z.^-1))./...
    ((1-p2*z.^-1).*(1-p1*z.^-1).*(1-p1'*z.^-1));
H{4}=@(z)-H{3}(z);
H{5}=@(z)c*((z.^-1-z2).*(z.^-1-z1).*(z.^-1-z1'))./...
    ((1-p2*z.^-1).*(1-p1*z.^-1).*(1-p1'*z.^-1));
H{6}=@(z)-H{5}(z);
H{7}=@(z)c*((1-z2*z.^-1).*(z.^-1-z1).*(z.^-1-z1'))./...
    ((1-p2*z.^-1).*(1-p1*z.^-1).*(1-p1'*z.^-1));
H{8}=@(z)-H{7}(z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% abs(H(exp(j*omega))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:8
    subplot(2,4,k)
    fplot(@(omega)abs(H{k}(exp(1j*omega))),[0,2*pi]);
    title(['|H_',num2str(k),'(e^{j*\omega})|'])
    xlabel('\omega')
end
%%%% The figures are identical,
% so we do not display them on the same graph.
%%%%%% zeros and poles
figure
zeros_H=[[z2;z1;z1'],[1/z2;z1;z1'],[1/z2;1/z1;...
    1/z1'],[z2;1/z1;1/z1']];
for k=1:4
    subplot(2,2,k)
    zplane(zeros_H(:,k),[p2;p1;p1']);
    title(['zeros and poles of ','H_',num2str(2*k-1),'(e^{j*\omega})'...
        ',H_',num2str(2*k),'(e^{j*\omega})'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
%%%%%%%%%%% phase response
for k=1:8
    subplot(2,4,k)
    fplot(@(omega)angle(H{k}(exp(1j*omega))),[0,2*pi]);
    title(['\angle H_',num2str(k),'(e^{j*\omega})'])
    xlabel('\omega')
end
%%%%%%%%%% group delay
denominator=conv([1,-p1],conv([1,-p1'],[1,-p2]));
H_nominator{1}=c*conv([1,-z2],conv([1,-z1'],[1,-z1]));
H_nominator{2}=-H_nominator{1};
H_nominator{3}=c*conv([-z2,1],conv([1,-z1'],[1,-z1]));
H_nominator{4}=-H_nominator{3};
H_nominator{5}=c*conv([-z2,1],conv([-z1',1],[-z1,1]));
H_nominator{6}=-H_nominator{5};
H_nominator{7}=c*conv([1,-z2],conv([-z1',1],[-z1,1]));
H_nominator{8}=-H_nominator{7};

figure
for k=1:8
    subplot(2,4,k)
    grpdelay(H_nominator{k},denominator);
    title(['\tau_{g',num2str(k),'}( H_',num2str(k),'(e^{j*\omega}))'])
    xlabel('\omega')
end


%%%%%%%%%% Pulse response
% Inverse Z transform of the denominator
a=1/(13/19-exp(-2*1j*pi/3)-6/19*exp(-1j*pi/3));
h=@(n)(n>=0).*(36/283*0.3.^n+2*real(a*p1.^n));
figure
m=0:50;
h_values=zeros(8,length(m));
for k=1:8
    subplot(2,4,k)
    h_values(k,:)=H_nominator{k}(1)*h(m)+H_nominator{k}(2)*h(m-1)+...
        H_nominator{k}(3)*h(m-2)+H_nominator{k}(4)*h(m-3);
    plot(m,h_values(k,:));
    title(['h_',num2str(k),'[n]'])
    xlabel('n')
end
%%%%%%%%%% cumulative energy
figure
cum_sum_h = cumsum(h_values.^2');
for k=1:8
    subplot(2,4,k)
    plot(m,cum_sum_h(:,k));
    title(['cumsum(h_',num2str(k),'[n])'])
    xlabel('n')
end
figure
%%%%%%%%%%% phase response
subplot(2,2,1)
hold all

for k=1:8
    fplot(@(omega)angle(H{k}(exp(1j*omega))),[0,2*pi]);
end
title('\angle H(e^{j*\omega})')
xlabel('\omega')
legend('H_1','H_2','H_3','H_4','H_5','H_6','H_7','H_8')
%%%%
subplot(2,2,2)
hold all

for k=1:8
    grpdelay(H_nominator{k},denominator);
end
title('\tau_{g}( H(e^{j*\omega}))')
xlabel('\omega')
legend('H_1','H_2','H_3','H_4','H_5','H_6','H_7','H_8')

%%%%
subplot(2,2,3)
hold all

for k=1:8
    plot(m,h_values(k,:));
    
end
title('h[n]')
xlabel('n')
legend('h_1','h_2','h_3','h_4','h_5','h_6','h_7','h_8')
%%%%
subplot(2,2,4)
hold all

for k=1:8
    plot(m,cum_sum_h(:,k));
end
title('cumsum(h[n])')
xlabel('n')
legend('h_1','h_2','h_3','h_4','h_5','h_6','h_7','h_8')


%%% Pay attention that H7,H8 are the minimal phase systems.
%They attain the lowest group delay and the maximal cummulative sum
