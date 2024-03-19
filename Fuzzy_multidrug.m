% New idea implementation for Multi-drug
% Fuzzy rule base + multi-drug equations from base paper (panjwani 2021)
% 
% Created by: Sajal Debnath (23th april, 2023)

%% 
a = 0.5;
b = 0.477;
c = 0.218;
d = 0.05;
e = 0.0173;
alpha = 0.1;
K = 1e9; % Initial Normal cell
gama_a = 0.32;
gama_b = 0.25;
miu = 0.4;

lambda_a = 0.0084; % irinotecan
lambda_b = 0.0092; % docetaxel

%% New parameter for multidrug
% ta = 2; %first dose time instance (irinotecan)
% tb = 3; %first dose time instance (docetaxel)
bita_a = 0.008; % Probability for mutation from sensitive population to resistant population for drug (irinotecan)
bita_b = 0.014; % Probability for mutation from sensitive population to resistant population for drug (docetaxel)


%%

P = zeros(50,50);    % Proliferating cells
Q = zeros(50,50);    % Quiescent cells
P_a = zeros(50,50);  % Proliferating cell resistance to drug (A-> irinotecan)
Q_a = zeros(50,50);  % Quiescent cell resistance to drug (A--> irinotecan)
P_b = zeros(50,50);  % resistance to docetaxel
Q_b = zeros(50,50);  % resistance to docetaxel
P_ab = zeros(50,50); % resistance of both drug
Q_ab = zeros(50,50); % resistance of both drug

M = zeros(50,50);
tem_M = zeros(8,8);
G_a = zeros(50,50);
G_b = zeros(50,50);
D_a = zeros(50,50);
D_b = zeros(50,50);
T = zeros(50,50); % toxicity is not needed for multidrug

%% 

P(1,1) = 2e11;
Q(1,1) = 8e11;
M(1,1) = 1e9; % Initial normal cell

Dose_a(1) = 0;
Dose_b(1) = 0;

D_a(1,1) = 0;
D_b(1,1) = 0;

map_rule = zeros(25,1);

T(1,1) = 0; % not needed

tem_D_a = zeros(1,8);
tem_D_b = zeros(1,8);
tem_G_a = zeros(1,8); 
tem_G_b = zeros(1,8);

%%

day_a = [2, 8, 15, 22, 29, 36, 43, 50, 57, 64, 71, 78, 85]; 
day_b = [2, 22, 43, 64, 85]; 

% Dose = [25, 7, 20, 20, 10, 20, 20, 10, 20, 10, 20, 20, 20];  % replaced with pattern7 and pattern 21
for i = 2:1:87
       
    input = [log10(P(1, i-1)) log10(M(1, i-1))];
    if(find(day_a == i))
%         [map_rule] = rule_map(map_rule, log10(P(1, i-1)), log10(M(1, i-1)));
        ut_a = evalfis(Pro_norm_fuzzy_16may, input);
    else
        ut_a = 0;
    end
    
    if(find(day_b == i))
%         [map_rule] = rule_map(map_rule, log10(P(1, i-1)), log10(M(1, i-1)));
        ut_b = evalfis(Pro_norm_fuzzy_16may, input);
    else
        ut_b = 0;
    end
    
%     total_dose = ut_a+ut_b;
%     if total_dose > 50
%         greater_part = total_dose-50;
%         have_to_less = greater_part/2;
%         ut_a = ut_a - have_to_less;
%         ut_b = ut_b - have_to_less;
%     end
    
    
    % if same day for both drug
    if ut_a>0 && ut_b>0
        tem_M(1,1) = M(1, i-1);
        tem_D_a(1,1) = D_a(1, i-1);
        tem_D_b(1,1) = D_b(1, i-1);
        tem_G_a(1,1) = 0; 
        tem_G_b(1,1) = 0; 
        UT_A = ut_a;
        UT_B = ut_b;
        for j=2:7
            
            tem_D_a(1, j) = tem_D_a(1, j-1) + UT_A - gama_a*tem_D_a(1, j-1);
            tem_D_b(1, j) = tem_D_b(1, j-1) + UT_B - gama_b*tem_D_b(1, j-1);
            
            tem_G_a(1, j) = lambda_a* tem_D_a(1,j);
            tem_G_b(1, j) = lambda_b* tem_D_b(1,j);
            tem_M(1, j) = tem_M(1, j-1) + alpha * tem_M(1, j-1) * (1- tem_M(1, j-1)/K) - (tem_G_a(1, j) + tem_G_b(1, j))*tem_M(1, j-1);
            
            UT_A = 0;
            UT_B = 0;
        end
        fprintf("%d) After 7 day Normal cell will be %.3e\n",i,tem_M(1, 7));
        if log10(tem_M(1, 7))<=8.14
            input_2 = [log10(tem_M(1, 7))];
            percent_dose_dec = evalfis(Normal_cell_to_percent_decrease_v1, input_2);
            fprintf("Dose should decrease: %.2f percent.\n",percent_dose_dec);
            
            ut_a = ut_a - ut_a*(percent_dose_dec/100);
            ut_b = ut_b - ut_b*(percent_dose_dec/100);
        end
    end
    
    
    Dose_a(i) = ut_a;
    Dose_b(i) = ut_b;
    
    D_a(1, i) = D_a(1, i-1) + ut_a - gama_a*D_a(1, i-1);
    D_b(1, i) = D_b(1, i-1) + ut_b - gama_b*D_b(1, i-1);
    
    % T(1, i) = T(1, i-1) + D(1, i) - miu*T(1, i-1);  % toxicity is not applied for multidrug
    G_a(1, i) = lambda_a* D_a(1,i);
    G_b(1, i) = lambda_b* D_b(1,i);
    
    T(1, i) = T(1, i-1) + D_a(1, i) + D_b(1, i) - miu*T(1, i-1); % toxicity
    
    P(1, i) = P(1, i-1) + (a-b-c)*(1- bita_a*ut_a -bita_b*ut_b)*P(1, i-1) + d*Q(1, i-1) - (G_a(1, i)+G_b(1, i))*P(1, i-1);          
    Q(1, i) = Q(1, i-1) + c*P(1, i) - (d+e)*Q(1, i-1);
    
    P_b(1, i) = P_b(1, i-1) + (a-b-c)*(1-bita_a*ut_a)*P_b(1, i-1) + d*Q_b(1, i-1) + bita_b*ut_b*P(1, i-1) - G_a(1, i)*P_b(1, i-1);          
    Q_b(1, i) = Q_b(1, i-1) + c*P_b(1, i) - (d+e)*Q_b(1, i-1);
    
    P_a(1, i) = P_a(1, i-1) + (a-b-c)*(1-bita_b*ut_b)*P_a(1, i-1) + d*Q_a(1, i-1) + bita_a*ut_a*P(1, i-1) - G_b(1, i)*P_a(1, i-1);          
    Q_a(1, i) = Q_a(1, i-1) + c*P_a(1, i) - (d+e)*Q_a(1, i-1);
    
    P_ab(1, i) = P_ab(1, i-1) + (a-b-c)*P_ab(1, i-1) + d*Q_ab(1, i-1) + bita_a*ut_a*P_b(1, i-1) + bita_b*ut_b*P_a(1, i-1);          
    Q_ab(1, i) = Q_ab(1, i-1) + c*P_ab(1, i) - (d+e)*Q_ab(1, i-1);
    
    M(1, i) = M(1, i-1) + alpha * M(1, i-1) * (1- M(1, i-1)/K) - (G_a(1, i) + G_b(1, i))*M(1, i-1);

end

Total_cancerous_cell = P(1,87) + Q(1,87)+ P_a(1,87) + Q_a(1,87) + P_b(1,87) + Q_b(1,87) + P_ab(1,87) + Q_ab(1,87);


%% Results

p_reduce = ((2e11 - P(1,87))*100)/2e11;
q_reduce = ((8e11 - Q(1,87))*100)/8e11;
formatSpec = 'Pro reduce %4.2f, Quies reduce %4.2f , Nor-> %.3e\n';
fprintf(formatSpec,p_reduce,q_reduce,M(1,87));


%% Figure Draw
figure
hold on
y =[1:size(M,2)];
plot(y,Dose_a(1,:), 'r-*');
plot(y,Dose_b(1,:), 'b-^');

D_a = D_a + D_b;
plot(y,D_a(1,:), 'k-');
% plot(y,D_b(1,:), 'g-');

% plot(y,Q,'g-v');
% plot(y,x4,'b-o');
% plot(y,x5,'k-diamond');

xlabel('Population No.','FontName','Times');
ylabel('Cell No.','FontName','Times');
legend({'Irinotecan','Docetaxel','Drug Concentration'},'FontName','Times','Location','best');
grid on
hold off

%% Another Figure 

figure
plot(y,T(1,:), 'g-');

xlabel('No. of Days','FontName','Times');
ylabel('Toxicity','FontName','Times');
% legend({'Irinotecan','Docetaxel','Drug Concentration'},'FontName','Times','Location','best');
grid on


% figure
% % fis = readfis('Pro_norm_fuzzy_16may');
% fis = readfis('Normal_cell_to_percent_decrease_v1');
% plotmf(fis,'output',1)
% xlim([-3, 40]);
% 
% grid on

%%

function [map_rule]=rule_map(map_rule, pro, normal)
    fprintf("pro--> %.2f, normal--> %.2f\n", pro, normal);
    
    if (pro >= 10.75 && pro <= 11.5) && (normal >= 8.75 && normal <= 9.25)
        map_rule(1) =  map_rule(1) + 1;
    end
    if (pro >= 10.75 && pro <= 11.5) && (normal >= 8.5 && normal <= 9)
        map_rule(2) =  map_rule(2) + 1;
    end
    if (pro >= 10.75 && pro <= 11.5) && (normal >= 8.25 && normal <= 8.75)
        map_rule(3) =  map_rule(3) + 1;
    end
    if (pro >= 10.75 && pro <= 11.5) && (normal >= 8 && normal <= 8.5)
        map_rule(4) =  map_rule(4) + 1;
    end
    if (pro >= 10.75 && pro <= 11.5) && (normal >= 7.75 && normal <= 8.25)
        map_rule(5) =  map_rule(5) + 1;
    end
    
    %
    
    if (pro >= 10.3 && pro <= 11) && (normal >= 8.75 && normal <= 9.25)
        map_rule(6) =  map_rule(6) + 1;
    end
    if (pro >= 10.3 && pro <= 11) && (normal >= 8.5 && normal <= 9)
        map_rule(7) =  map_rule(7) + 1;
    end
    if (pro >= 10.3 && pro <= 11) && (normal >= 8.25 && normal <= 8.75)
        map_rule(8) =  map_rule(8) + 1;
    end
    if (pro >= 10.3 && pro <= 11) && (normal >= 8 && normal <= 8.5)
        map_rule(9) =  map_rule(9) + 1;
    end
    if (pro >= 10.3 && pro <= 11) && (normal >= 7.75 && normal <= 8.25)
        map_rule(10) =  map_rule(10) + 1;
    end
    
    % 
    
    if (pro >= 10 && pro <= 10.75) && (normal >= 8.75 && normal <= 9.25)
        map_rule(11) =  map_rule(11) + 1;
    end
    if (pro >= 10 && pro <= 10.75) && (normal >= 8.5 && normal <= 9)
        map_rule(12) =  map_rule(12) + 1;
    end
    if (pro >= 10 && pro <= 10.75) && (normal >= 8.25 && normal <= 8.75)
        map_rule(13) =  map_rule(13) + 1;
    end
    if (pro >= 10 && pro <= 10.75) && (normal >= 8 && normal <= 8.5)
        map_rule(14) =  map_rule(14) + 1;
    end
    if (pro >= 10 && pro <= 10.75) && (normal >= 7.75 && normal <= 8.25)
        map_rule(15) =  map_rule(15) + 1;
    end
    
    %
    
    if (pro >= 9 && pro <= 10.3) && (normal >= 8.75 && normal <= 9.25)
        map_rule(16) =  map_rule(16) + 1;
    end
    if (pro >= 9 && pro <= 10.3) && (normal >= 8.5 && normal <= 9)
        map_rule(17) =  map_rule(17) + 1;
    end
    if (pro >= 9 && pro <= 10.3) && (normal >= 8.25 && normal <= 8.75)
        map_rule(18) =  map_rule(18) + 1;
    end
    if (pro >= 9 && pro <= 10.3) && (normal >= 8 && normal <= 8.5)
        map_rule(19) =  map_rule(19) + 1;
    end
    if (pro >= 9 && pro <= 10.3) && (normal >= 7.75 && normal <= 8.25)
        map_rule(20) =  map_rule(20) + 1;
    end
    
    %
    
    if (pro >= 9 && pro <= 10) && (normal >= 8.75 && normal <= 9.25)
        map_rule(21) =  map_rule(21) + 1;
    end
    if (pro >= 9 && pro <= 10) && (normal >= 8.5 && normal <= 9)
        map_rule(22) =  map_rule(22) + 1;
    end
    if (pro >= 9 && pro <= 10) && (normal >= 8.25 && normal <= 8.75)
        map_rule(23) =  map_rule(23) + 1;
    end
    if (pro >= 9 && pro <= 10) && (normal >= 8 && normal <= 8.5)
        map_rule(24) =  map_rule(24) + 1;
    end
    if (pro >= 9 && pro <= 10) && (normal >= 7.75 && normal <= 8.25)
        map_rule(25) =  map_rule(25) + 1;
    end
    
    
   
    
    
    
    if (pro >= 9 && pro <= 10)
        fprintf("Proliferating --> Very Small.\n");
    end
    if (pro >= 9 && pro <= 10.3)
        fprintf("Proliferating --> Small.\n");
    end
    if (pro >= 10 && pro <= 10.75)
        fprintf("Proliferating --> Medium.\n");
    end
    if (pro >= 10.3 && pro <= 11)
        fprintf("Proliferating --> Big.\n");
    end
    if (pro >= 10.75 && pro <= 11.5)
        fprintf("Proliferating --> Very Big.\n");
    end
    
    
%     check Normal cells
    
    
    if (normal >= 7.75 && normal <= 8.25)
        fprintf("Normal --> Very Small.\n");
    end
    if (normal >= 8 && normal <= 8.5)
        fprintf("Normal --> Small.\n");
    end
    if (normal >= 8.25 && normal <= 8.75)
        fprintf("Normal --> Medium.\n");
    end
    if (normal >= 8.5 && normal <= 9)
        fprintf("Normal --> Big.\n");
    end
    if (normal >= 8.75 && normal <= 9.25)
        fprintf("Normal --> Very Big.\n");
    end
    
    fprintf("\n");
end

















