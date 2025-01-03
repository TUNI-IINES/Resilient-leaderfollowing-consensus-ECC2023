

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%Resilient Leader following consensus%%%%%%%
%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%

clear all;
close all;
clc;

n = 5;
x_ini = [1 1.5 2 -1 3]';
% x_ini(5,1) = 1;
z_ini =  [1 1.7 2 -1 3]';
% z_ini(1,1) = 1;
d_ini = zeros(n,1);
d2_ini = zeros(n,1);
chi_ini =[x_ini ; z_ini ; d_ini; d2_ini ];
tspan = [0 10];
options = odeset('RelTol',1e-6);

[t chi] = ode45(@x_position, tspan, chi_ini, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta1 = 40;
beta2 = 10;
alpha = 0.1*((sin( t) + cos(2*t) + sin(3*t)))+1;
% I2_i1 = chi(:,3)+ chi(:, 13);
% I2_i2 = beta1*(alpha .* chi(:,8)) + chi(:, 18);
% I2_i3 = (beta1*(alpha .* chi(:,3)))+(beta2 .* chi(:,8)) + chi(:, 13) + chi(:, 18);
% gain1 = 1 ./(beta1 .* alpha);
% gain2 = beta2 .* gain1;
% I2_i1_est = gain1 .* (I2_i3 - (gain2 .* I2_i2));
% error = I2_i1 - I2_i1_est;
% % 
% I3_i1 = chi(:,4);
% I3_i2 = beta1*(alpha .* chi(:,8)) ;
% I3_i3 = (beta1*(alpha .* chi(:,4)))+(beta2 .* chi(:,8)) ;
% gain1 = 1 ./(beta1 .* alpha);
% gain2 = beta2 .* gain1;
% I3_i1_est = gain1 .* (I3_i3 - (gain2 .* I3_i2));
% error = I3_i1 - I3_i1_est;
% % % % 
% 
% 
% 
% sqr_err = (error).^2;

% %%%%%%%%%%%%%%%%%% Privacy Preservation %%%%%%%%%%%%%%%%%%%
I = chi(:,1);
I0_1a = chi(:,1) - beta1 * (alpha .* chi(:,6));
I0_1b = (beta2 * chi(:,6)) + (beta1 * (alpha .* chi(:,1)));

%%%%%%%%%%%%%%%%%% Privacy Preservation %%%%%%%%%%%%%%%%%%%
% I = chi(:,4);
% I3_1a = chi(:,4) - beta1 * (alpha .* chi(:,9));
% I3_1b = (beta2 * chi(:,9)) + (beta1 * (alpha .* chi(:,4)));

n = 5;
s1 = [1 1 2 3 3 4 5];
s2 = [2 3 1 1 4 3 1];
G = digraph(s1,s2);
opts.Colors     = get(groot,'defaultAxesColorOrder');
opts.saveFolder = 'img/';
opts.width      = 17;
opts.height     = 5;
opts.fontType   = 'Times';
opts.fontSize   = 9;

% create new figure
fig = figure; clf

% scaling
fig.Units               = 'centimeters';
fig.Position(3)         = opts.width;
fig.Position(4)         = opts.height;
axis tight
% remove unnecessary white space
set(gca,'LooseInset',max(get(gca,'TightInset'), 0.02))
% export to png
fig.PaperPositionMode   = 'auto';
% plot(t, chi(:,1:5),'linewidth',1.5);
% ylabel('$x_i$ for all $i \in \mathcal V$', 'Interpreter','latex','FontSize',16);

% plot(t, abs(error),'linewidth',1.5)
% xlim([0 10])
% ylim([-1 1])
% hold on 
% plot(t, I3_i1,'linewidth',1.5)

% %%%%%%%%%%%%%Privacy %%%%%%%%%%%%%%%%%%%%%
plot(t, I,'linewidth',1.5);
hold on
plot(t, I0_1a,'linewidth',1.5);
plot(t, I0_1b,'linewidth',1.5);
legend('$x_0$', '$I_{0,a}$','$I_{0,b}$','Interpreter','latex','FontSize',16)
%%%%%%%%%%%%%Privacy %%%%%%%%%%%%%%%%%%%%%
% plot(t, I,'linewidth',1.5);
% plot(t, I3_1a,'linewidth',1.5);
% hold on
% 
% plot(t, I3_1b,'linewidth',1.5);
%legend( '$I_{3,a}$','$I_{3,b}$','Interpreter','latex','FontSize',16)
% % H = plot(G,'EdgeColor','r', 'linewidth',1.5);
% % highlight(H,[1 2 3 4],'NodeColor','r')
% % color = get(fig,'Color');
% % set(gca,'XColor',color,'YColor',color,'TickDir','out')
%%%%%%%%%%% AttackedLink%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend('$|{I}^2_{1,1} - \hat{I}^2_{1,1}|$', '${I}^2_{1}$', 'Interpreter','latex','FontSize',16)
% legend('$|\tilde{I}^3_{1,1} - \hat{I}^3_{1,1}|$', 'Interpreter','latex')

%%%%%%%%%%% No Attacked Link%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ylabel('$|{I}^2_{1,1} - \hat{I}^2_{1,1}|$', 'Interpreter','latex','FontSize',16);
% ylabel('$|{I}^3_{1,1} - \hat{I}^3_{1,1}|$', 'Interpreter','latex','FontSize',16);
xlabel('Time (Seconds)','FontSize',16);
% ylabel('$x_3 $' , 'Interpreter','latex','FontSize',16);
% save2pdf('privacy2_low_beta2',fig,600);
save2pdf('privacy1_high_beta2',fig,600);

function chi_dot = x_position(t, chi)
n = 5;
p = 10;
s1 = [1 2 2 3 4 4 5];
s2 = [2 3 4 2 2 5 4];
G = digraph(s1,s2);


in_deg = [];
for kk = 1 : n
    dd = find(s2 == kk);
    card_dd = length(dd);
    in_deg = [in_deg card_dd];
end
Adj = adjacency(G);
D = diag(in_deg);
L = D - Adj'; 


F1 = -eye(n,n);
F2 = -2*eye(n,n);
Ba = [0 0 0 0 0;
      0 2 -1 0 0;
      0 -1 0 0 0;
      0 0 0 0 0;
      0 0 0 0 1];

beta1 = 40;
beta2 = 10;
  alpha = 0.1*((sin( t) + cos(2*t) + sin(3*t)))+1;
% alpha = 1;
for i = 1 : n
%     if i == 1 
%         g1 = 0;
%     else
%         g1 = 1;
%     end
%     x_dot(i,1) =  - (L(i,:) * chi(1:5,1)) + (alpha * beta1 * (L(i,:) * chi(6:10,1)))  + g1*chi((2*n)+i,1) ;
%     z_dot(i,1) = - (alpha * beta1 * (L(i,:) * chi(1:5,1))) - (beta2 * L(i,:) * chi(6:10,1)) + g1*chi(i+(3*n),1);

%     x_dot(i,1) =  - (L(i,:) * chi(1:5,1)) + (alpha * beta1 * (L(i,:) * chi(6:10,1)))  + chi((2*n)+i,1) ;
%     z_dot(i,1) = - (alpha * beta1 * (L(i,:) * chi(1:5,1))) - (beta2 * L(i,:) * chi(6:10,1)) + chi(i+(3*n),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ForPrivacy%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 x_dot(i,1) =  - (L(i,:) * chi(1:5,1)) + (alpha * beta1 * (L(i,:) * chi(6:10,1)))  + chi((2*n)+i,1);
 z_dot(i,1) = - (alpha * beta1 * (L(i,:) * chi(1:5,1))) - (beta2 * L(i,:) * chi(6:10,1)) + chi(i+(3*n),1);

end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d_dot = F1 * chi(11:15,:) + Ba * chi(1:5,:) ;
d2_dot = F2 * chi(16:20,:) + Ba * chi(1:5,:) ;
chi_dot = [x_dot ; z_dot; d_dot; d2_dot];

end