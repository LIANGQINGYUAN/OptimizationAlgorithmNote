function [winner,bestfitness] = gar(Domain,N,fit_func,options)
% function winner = GAR(Domain,N,fit_func)
% Function call: GAR(Domain,N,¡¯f¡¯)
% Domain = search space; e.g., [-2,2;-3,3] for the space [-2,2]x[-3,3]
% (number of rows of Domain = dimension of search space)
% N = population size (must be an even number)
% f = name of fitness value function
%
%Options:
%print = options(1);
%selection = options(5);
%max_iter=options(14);
%p_c = options(18);
%p_m = p_c/100;
%
%Selection:
% options(5) = 0 for roulette wheel, 1 for tournament

clf;
if nargin ~= 4
options = [];
if nargin ~= 3
disp('Wrong number of arguments.');
return;
end
end
if length(options) >= 14
if options(14)==0
options(14)=3*N;
end
else
options(14)=3*N;
end
if length(options) < 18
options(18)=0.75; %optional crossover rate
end
%format compact;
%format short e;
options = foptions(options);
print = options(1);
selection = options(5);
max_iter=options(14);
p_c = options(18);
p_m = p_c/100;
n = size(Domain,1);
lowb = Domain(:,1)';
upb = Domain(:,2)';
bestvaluesofar = 0;
for i = 1:N,
P(i,:) = lowb + rand(1,n).*(upb-lowb);
%Initial evaluation
fitness(i) = feval(fit_func,P(i,:));
end
[bestvalue,best] = max(fitness);
if bestvalue > bestvaluesofar,
bestsofar = P(best,:);
bestvaluesofar = bestvalue;
end

for k = 1:max_iter,
%Selection
fitness = fitness - min(fitness); % to keep the fitness positive
if sum(fitness) == 0,
disp('Population has identical chromosomes -- STOP');
disp('Number of iterations:');
disp(k);
for i = k:max_iter,
myupper(i)=myupper(i-1);
average(i)=average(i-1);
mylower(i)=mylower(i-1);
end
break;
else
fitness = fitness/sum(fitness);
end
if selection == 0,
%roulette-wheel
cum_fitness = cumsum(fitness);
for i = 1:N,
tmp = find(cum_fitness-rand>0);
m(i) = tmp(1);
end
else
%tournament
for i = 1:N,
fighter1=ceil(rand*N);
fighter2=ceil(rand*N);
if fitness(fighter1)>fitness(fighter2),
m(i) = fighter1;
else
m(i) = fighter2;
end
end
end
M = zeros(N,n);
for i = 1:N,
M(i,:) = P(m(i),:);
end

%Crossover
Mnew = M;
for i = 1:N/2
ind1 = ceil(rand*N);
ind2 = ceil(rand*N);
parent1 = M(ind1,:);
parent2 = M(ind2,:);
if rand < p_c
a = rand;
offspring1 = a*parent1+(1-a)*parent2+(rand(1,n)-0.5).*(upb-lowb)/10;
offspring2 = a*parent2+(1-a)*parent1+(rand(1,n)-0.5).*(upb-lowb)/10;
%do projection
for j = 1:n,
if offspring1(j)<lowb(j),
offspring1(j)=lowb(j);
elseif offspring1(j)>upb(j),
offspring1(j)=upb(j);
end
if offspring2(j)<lowb(j),
offspring2(j)=lowb(j);
elseif offspring2(j)>upb(j),
offspring2(j)=upb(j);
end
end
Mnew(ind1,:) = offspring1;
Mnew(ind2,:) = offspring2;
end
end
%Mutation
for i = 1:N,
if rand < p_m,
a = rand;
Mnew(i,:) = a*Mnew(i,:) + (1-a)*(lowb + rand(1,n).*(upb-lowb));
end
end
P = Mnew;
%Evaluation
for i = 1:N,
fitness(i) = feval(fit_func,P(i,:));
end
[bestvalue,best] = max(fitness);
if bestvalue > bestvaluesofar,
bestsofar = P(best,:);
bestvaluesofar = bestvalue;
end
myupper(k) = bestvalue;
average(k) = mean(fitness);
mylower(k) = min(fitness);
end %for
if k == max_iter,
disp('Algorithm terminated after maximum number of iterations:');
disp(max_iter);
end
winner = bestsofar;
bestfitness = bestvaluesofar;
if print,
iter = [1:max_iter]';
plot(iter,myupper,'o:',iter,average,'x-',iter,mylower,'*--');
legend('Best', 'Average', 'Worst');
xlabel('Generations','Fontsize',14);
ylabel('Objective Function Value','Fontsize',14);
set(gca,'Fontsize',14);
hold off;
end

function y=f_wave(x);
y=x(1)*sin(x(1)) + x(2)*sin(5*x(2));
