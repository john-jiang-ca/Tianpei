function [ offspring1, offspring2 ] = crossWalk( parent1, parent2 )
%lthis subroutine implement crosswalk to generate two offsprings by two
%parents
Nt=length(parent1);
offspring1=parent1;
offspring2=parent2;
R1=randi(Nt,1);
offspring_tmp=offspring1([1:R1]);
offspring1([1:R1])=offspring2([1:R1]);
offspring2([1:R1])=offspring_tmp;

end

