


figure
subplot(211)
plot(Xtest(1:10,:)')


subplot(212)
plot(mean(abs(theta.A(:,1,:)),3) / mean(mean(abs(theta.A(:,1,:)),3)))
hold on
plot(mean(abs(theta.A(:,2,:)),3) / mean(mean(abs(theta.A(:,2,:)),3)))
plot(mean(abs(theta.A(:,3,:)),3) / mean(mean(abs(theta.A(:,3,:)),3)))
 
% subplot(122)
% plot(std(squeeze(theta.A(:,1,:))') / mean(std(squeeze(theta.A(:,1,:))')))
% hold on
% plot(std(squeeze(theta.A(:,2,:))') / mean(std(squeeze(theta.A(:,2,:))')))
% plot(std(squeeze(theta.A(:,3,:))') / mean(std(squeeze(theta.A(:,3,:))')))

