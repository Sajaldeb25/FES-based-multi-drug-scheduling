%{
    --> Appling PSO multiobjective algorithm with Panetta(1995) equations.
    --> Sajal Debnath (Feb 2023)
%}


%%
close all
clc
format long
warning off
global No Co taug rowg keff gamaa Tmax Ccum Dmax tf interv v xp xg dose flc Cmax Cth
interv= 14;
taug=150;
rowg=1e12;
No=1e10;
Co=0;
keff=2.7e-2;
Ccum=4.1e3;
tf=364;
gamaa=0.27;
c=2;
w=0.7;
particles=50;
iteration=50;
var=13;
c_cf=0;
changed_dose = zeros(50,13);

%initialization


varlo=[15 15 15 15 15 15 15 15 15 15 15 15 15];
varhi=[23 23 23 23 23 23 23 23 23 23 23 23 23];


for m=1:particles
    for n=1:var
        v(m,n)=0; %velocity particles
        a=varlo(n);
        b=varhi(n);
        %disp(a);
        
        x(m,n)=a+rand*(b-a); %position particles
        xp(m,n)=x(m,n);
    end
end
% disp(xp);
for ixno = 1:particles
  [cost(ixno,:)] = obj_fitness_value(xp(ixno,:)); %fitness value
end



%% find the best value
   
[best_performance,location]=min( cost);
fg=best_performance;
xg=x(location,: );
for gen=1:iteration
    for m=1:particles
        for n=1:var
            v(m,n)=(w*v(m,n))+(c*rand*(xp(m,n)-x(m,n)))+(c*rand*(xg(n)-x(m,n)));
            x(m,n)=x(m,n)+v(m,n);
        end
%         disp(x);
    end
    
    for m=1:particles
        for i=1:var
          if x(m,i)< varlo(i)
             x(m,i)= varlo(i);
          end
          if x(m,i)> varhi(i)
             x(m,i)= varhi(i);
          end
        end
    end
        
    for ixno = 1:particles
        [cost1(ixno,:), poli(ixno), quo(ixno), nor(ixno), changed_dose(ixno,:)] = obj_fitness_value(x(ixno,:)); %fitness value
    end  
   
    x = changed_dose;
    %% print optimal result
    [value, pos] = min(poli);
    p_reduce = ((2e11 - value)*100)/2e11;
    q_reduce = ((8e11 - quo(pos))*100)/8e11;
    nor_cell = nor(pos);
    formatSpec = 'Gen->%d, pro reduce %4.2f, Quies reduce %4.2f , Nor-> %.3e\n';
    fprintf(formatSpec,gen,p_reduce,q_reduce,nor_cell);
    Global_dose(gen,:) = x(pos,:);
    
    %%
    for ixno = 1:particles
        if cost1(ixno,:)<cost(ixno,:)
             cost(ixno,:)=cost1(ixno,:);
             xp(m,n)=x(m,n);
        else
            xp(m,n)=xp(m,n);
        end
    end
    [B_fg,location]=min(cost);
    
    %COMPARE global
    if B_fg<fg
        fg=B_fg;
        xg=xp(location,:);
    end
    
    c_cf=c_cf+1;
    best_cf_ac=fg;
    
end
   
min_cost=fg;

t_cf=1:c_cf;
%figure
figure(1)
plot(cost(:,1),cost(:,2),'b*');
title(['Iterations--', num2str(iteration)]);
xlabel('Z1');ylabel('Z2');

shg
drawnow;

%%
% Gen->44, pro reduce 92.34, Quies reduce 89.15 , Nor-> 1.842e+08
% dose = [23,20.57,22.28,15,20.52,23,15.41,16.67,16.23,15.33,20.61,18.04,23]





