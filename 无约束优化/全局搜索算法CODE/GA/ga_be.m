function [winner,bestfitness] = ga_be(L,N,fit_func,options)
% function winner = GA(L,N,fit_func)
% Function call: GA(L,N,¡¯f¡¯)
% L = length of chromosomes
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

P = rand(N,L)>0.5;
bestvaluesofar = 0;

%Initial evaluation
for i = 1:N,
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
    M = zeros(N,L);
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
            crossover_pt = ceil(rand*(L-1));
            offspring1 = [parent1(1:crossover_pt) parent2(crossover_pt+1:L)];
            offspring2 = [parent2(1:crossover_pt) parent1(crossover_pt+1:L)];
            Mnew(ind1,:) = offspring1;
            Mnew(ind2,:) = offspring2;
        end
    end

    %Mutation
    mutation_points = rand(N,L) < p_m;
    P = xor(Mnew,mutation_points);
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

function y=f_manymax(x);
y=-15*(sin(2*x))^2-(x-2)^2+160;
function y=fit_func1(binchrom);
%1-D fitness function
f='f_manymax';
range=[-10,10];
x=bin2dec(binchrom,range);
y=feval(f,x);

function y=f_peaks(x);
y=3*(1-x(1)).^2.*exp(-(x(1).^2)-(x(2)+1).^2) -10.*(x(1)/5-x(1).^3-x(2).^5).*exp(-x(1).^2-x(2).^2) - exp(-(x(1)+1).^2-x(2).^2)/3;
function y=fit_func2(binchrom);
%2-D fitness function
f='f_peaks';
xrange=[-3,3];
yrange=[-3,3];
L=length(binchrom);
x1=bin2dec(binchrom(1:L/2),xrange);
x2=bin2dec(binchrom(L/2+1:L),yrange);
y=feval(f,[x1,x2]);