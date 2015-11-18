function [ symOut ] = RSVD_GA( y, H, symGeneration_prev, M, symConstell, C, tol, epsilon, Pm, S_AREA, Alpha, SNRd, pav, population, generation, maxGeneration )
%hybrid algorithm of support vector detector and Genetic algorithm
% 
Np=population;   %size of the population
Ng=maxGeneration;   %size of the generation
Nc=Np/2;  %number of the chromosomes to do crosswalk 
pc=0.95; %crosswalk probability
pm=Pm; %mutation probability
alpha=Alpha; %size of the population to do PIC
S_area=S_AREA; %searching area of the initial population
Nr=length(H(:,1)); %number of receive antennas
Nt=length(H(1,:));  %number of transmit antennas
noiseV=1/SNRd;  %the variance of noise
Threshold=Nr*noiseV+2*sqrt(Nr*noiseV^(2)); %the threshold according to non-central Chi-square distribution
%% Generate initial generation of chromosomes based on the preliminary estimation from RSVD
if (generation==0)
H_r=[real(H), -imag(H); imag(H), real(H)];
y_r=[real(y);imag(y)];
[ symOut_current] = RSVR( H_r,  y_r, SNRd,  M, pav,C,tol, epsilon );
symGeneration=zeros(Nt,Np);
symGeneration(:,1)=symOut_current;
for count=2:Np   %generate initial population by perturbing the output of RSVD
    symGeneration_tmp=symOut_current;
%     R1=randi(Nt,1);
    R1=randperm(Nt,S_area);
    for count1=1:S_area
    R2=randi(M,1);
    while(symGeneration_tmp(R1(count1))==symConstell(R2))
        R2=randi(M,1);
    end
    symGeneration_tmp(R1(count1))=symConstell(R2);
    end
    symGeneration(:,count)=symGeneration_tmp;
end
symGeneration_prev=symGeneration;
% symOut_prev=symOut_current;
end
symGeneration_current=symGeneration_prev;

%% Genetic Selection
fitness=zeros(Np,1);
fitness_matrix=y*ones(1,Np)-H*symGeneration_current;
for count=1:Np
    fitness(count)=norm(fitness_matrix(:,count))^(2);  %calculate the fitness function 
end
[value, order]=sort(fitness, 'ascend');
if(value(1)<Threshold)
    symOut=symGeneration_current(:,order(1));
    return
end
% parent1=symGeneration_current(:,order(1));  %elite are chosen as the parents for crosswalk
% parent2=symGeneration_current(:,order(2));
% symGeneration_current(:,Np-1)=zeros(Nt,1);   %the weakest ones are chosen to be sacrified
% symGeneration_current(:,Np)=zeros(Nt,1);
%% Crosswalk
% if (rand(1)<=pc)   %obey the probability to perform crosswalk
%    R1=randi(Nt,1);
%    parent_tmp=parent1([1:R1]);
%    parent1([1:R1])=parent2([1:R1]);
%    parent2([1:R1])=parent_tmp;
% end
first=1;
last=Nc;
while (first<last)
    [symGeneration_current(:,order(first+Nc)), symGeneration_current(:,order(last+Nc))]=...
        crossWalk(symGeneration_current(:,order(first)), symGeneration_current(:,order(last)));
    first=first+1;
    last=last-1;
end



%% Mutation
for count=1:Nc
    R1=rand(1);
    offspring_tmp=symGeneration_current(:,order(count+Nc));
    if(R1<=pm)
        R2=randi(Nt,1);
         R3=randi(M,1);
        while(offspring_tmp(R2)==symConstell(R3))
            R3=randi(M,1);
        end
        offspring_tmp(R2)=symConstell(R3);     
    end
end

%% Parallel Interference Cancellation
Npic=floor(alpha*Np);
list=randperm(Np);
stage=1;
maxStage=4;
G=H'*H;
[value1, order1]=sort(diag(G), 'ascend');
W=eye(Nt);
for count=1:Npic
% symOut_prev=symGeneration_current(:,order(count));    %choose the worst ones in the population to do local search 
%    symGeneration_current(:, order(count))=PIC_GA(y, H, symOut_prev, SNRd, M, pav, stage, maxStage, symConstell, order1,W);
symOut_prev=symGeneration_current(:,list(count));    %choose the worst ones in the population to do local search 
   symGeneration_current(:, list(count))=PIC_GA(y, H, symOut_prev, SNRd, M, pav, stage, maxStage, symConstell, order1, W);
end
%% stopping when the number of generation reaches the maximum
if (generation==maxGeneration)
    fitness=zeros(Np,1);
    fitness_matrix=y*ones(1,Np)-H*symGeneration_current;
for count=1:Np
    fitness(count)=norm(fitness_matrix(:,count))^(2);  %calculate the fitness function 
end
[value, index]=min(fitness);
    symOut=symGeneration_current(:,index);
    return 
end
%% Next generation 
 [symOut]=RSVD_GA( y, H, symGeneration_current, M, symConstell, C, tol, epsilon, pm, S_area, alpha, SNRd,pav, population, generation+1, maxGeneration );
 
end

