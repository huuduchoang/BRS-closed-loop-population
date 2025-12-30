x=-10:.001:10;
figure
plot(x,max(0,x),'LineWidth',3)
axis([-10 10 -.5 1])
xlabel('x')
ylabel('ReLu(x)')
set(gca,'FontSize',20)
grid on
hold on
softplus=@(x,r)log(1+exp(x*r))/r;
for r=1:2:40
%plot(x,(x+1/r).*exp(r*x)./(1+exp(r*x)),'LineWidth',2)
plot(x,softplus(x,r),'LineWidth',2)
end
