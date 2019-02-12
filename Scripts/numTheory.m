s = 0:0.01:1.5;
ap = 0.5*(1-i*s/2 - s.^2) + 0.5*sqrt(1-i*3*s-9/4*s.^2 + i*s.^3 + s.^4);
am = 0.5*(1-i*s/2 - s.^2) - 0.5*sqrt(1-i*3*s-9/4*s.^2 + i*s.^3 + s.^4);

figure(1);
hold off;
plot(s,abs(ap),'k');
hold on;
plot(s,abs(am),'k--');
plot(s,ones(size(s)),'k:');
plot(sqrt(2)*ones(size(s)),s,'k:');
hold off;
axis([0 1.5 0 1.4+eps])
title('|a_\pm| vs. \omega\Deltat');
xlabel('\omega\Deltat');
ylabel('|a_\pm|','Rotation',pi/2);
legend('|a_+|','|a_-|','Location','NorthWest');
text('Interpreter','latex',...
       'String','$$\sqrt{2} \rightarrow$$',...
       'HorizontalAlignment','right',...
       'Position',[sqrt(2) 0.2],...
       'FontSize',12);

amAngle = angle(am);
amAngle(1) = pi/2;
       
figure(2);
hold off;
plot(s,angle(ap),'k');
hold on;
plot(s(1:142),amAngle(1:142),'k--');
plot(s(143:151),amAngle(143:151),'k--');
plot(s,zeros(size(s)),'k:');
plot(sqrt(2)*ones(size(-4:0.1:4)),(-4:0.1:4),'k:');
hold off;
axis([0 max(s) -pi-eps pi+eps]);
title('arg(a_\pm) vs. \omega\Deltat');
xlabel('\omega\Deltat');
ylabel('arg(a_\pm)');
legend('|a_+|','|a_-|','Location','NorthWest');
text('Interpreter','latex',...
       'String','$$\sqrt{2} \rightarrow$$',...
       'HorizontalAlignment','right',...
       'Position',[sqrt(2) 1],...
       'FontSize',12);
