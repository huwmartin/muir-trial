clear;
close all;

%%% Parameters
CrossoverProbability = 1;
MutationProbability = 0.8;
Generations = 300;
PopulationSize = 30;
SequenceLength = 30;
Environment = 'muir_world.txt';

Population = zeros(PopulationSize, SequenceLength);

%%% Generate random population
for i = 1:PopulationSize
    TempChromosome = zeros(1, SequenceLength);
    
    for j = 1:SequenceLength/3
        TempChromosome((j*3)-2) = randi([1 4]);
        TempChromosome((j*3)-1) = randi([0 9]);
        TempChromosome(j*3) = randi([0 9]);
    end
    
    Population(i, :) = TempChromosome;
end

%%% Include extra column at end for scores
Population = [Population zeros(PopulationSize, 1)];

BestFitnessHistory = zeros(Generations, 1);
%%%     FitnessHistory = zeros(Generations, 3);

%%% Repeat `Generations` times to create a new population
for i = 1:Generations
    %%% Evaluate fitness scores
    for j = 1:PopulationSize
        [Fitness, Path] = simulate_ant(Environment, Population(j, 1:SequenceLength));
        Population(j, SequenceLength+1) = Fitness;
    end
    
    %%% Sort chromosomes in increasing fitness
    Population = sortrows(Population, SequenceLength+1);
    
    %%% Record best fitness in the population
    MaxFitness = max(Population(:, end));
    
    BestFitnessHistory(i) = MaxFitness;
    
    %%% Uncomment to record mean and weakest fitness as well 
    %%% as best (along with lines 31 and 73-77) 
    %%%     FitnessHistory(i, 1) = MaxFitness;
    %%%     FitnessHistory(i, 2) = mean(Population(:, end));
    %%%     FitnessHistory(i, 3) = min(Population(:, end));
    
    %%% Uncomment to use rank selection
    %%%     %%% Rank population instead of using fitness
    %%%     for j = 1:PopulationSize
    %%%         Population(j, SequenceLength+1) = j;
    %%%     end

    %%% Normalise fitness and calculate cumulative fitness variables to facilitate the selection step
    NormalisedFitness = Population(:, end)/sum(Population(:, end));
    CumulativeNormalisedFitness = cumsum(NormalisedFitness);
    
    if mod(i, 1) == 0
        fprintf('Iter. %d: Best solution: Fitness = %f\n', i, MaxFitness);
        figure(1);
        plot(BestFitnessHistory(1:i));
        xlabel('Iteration');
        ylabel('Fitness (best)');
        
        %%% Uncomment to plot best, mean, worst
        %%%     figure(2);
        %%%     plot(FitnessHistory(1:i, :));
        %%%     legend('Best', 'Mean', 'Worst')
        %%%     xlabel('Iteration');
        %%%     ylabel('Fitness (best)');
        pause(0.01);
    end
   
    %%% Elite, keep best 2
    PopulationNew = zeros(PopulationSize, SequenceLength);
    PopulationNew(1:2, :) = Population(PopulationSize-1:PopulationSize, 1:SequenceLength);
    PopulationNewSize = 2;
    
    %%% Repeat until new population is full
    while (PopulationNewSize < PopulationSize)
        %%% Use roulette wheel selection selection and pick two chromosomes
        %%% `rand` generates a (pseudo-)random number drawn from
        %%% the uniform distribution on the interval (0,1)
        R1 = rand;
        PIdx1 = find(CumulativeNormalisedFitness>=R1, 1, 'first');
        TempChromosome_1 = Population(PIdx1, 1:SequenceLength);
        
        R2 = rand;
        PIdx2 = find(CumulativeNormalisedFitness>=R2, 1, 'first');
        TempChromosome_2 = Population(PIdx2, 1:SequenceLength);
        
        %%% Perform crossover of `TempChromosome_1` and `TempChromosome_2` with probability, `CrossoverProbability`:
        if (rand < CrossoverProbability)
            %%% `SequenceLength/3` ensures that the cross point falls at the start of an instruction set
            %%% We `-2` to ensure that the cross point can at most be at the 8th instruction set so that 
            %%% the second cross point can be at the end of this set (so that there is a two-point crossover)
            CrossPoint_1 = randi([1 ((SequenceLength/3)-2)])*3;
            CrossPoint_2 = randi([(CrossPoint_1/3) ((SequenceLength/3)-1)])*3;

            TempChromosome_3 = [TempChromosome_1(1:CrossPoint_1), TempChromosome_2((CrossPoint_1+1):CrossPoint_2), TempChromosome_1(CrossPoint_2+1:SequenceLength)];
            TempChromosome_4 = [TempChromosome_2(1:CrossPoint_1), TempChromosome_1((CrossPoint_1+1):CrossPoint_2), TempChromosome_2(CrossPoint_2+1:SequenceLength)];

            TempChromosome_1 = TempChromosome_3;
            TempChromosome_2 = TempChromosome_4;
        end
        
        %%% Perform mutation on `TempChromosome_1` and `TempChromosome_2` with probability, `MutationProbability`:
        %%% For each chromosome, randomly picking two mutation points (but same instruction digit) and 
        %%% swap the corresponding values
        if (rand < MutationProbability)
            %%% Mutation by assigning a new random value
            TempChromosome_1 = random_value_mutation(SequenceLength, TempChromosome_1);
            
            %%% Mutation by order changing (also uncomment line 129 and comment out lines 119 and 126) 
            %%%     TempChromosome_1 = order_changing_mutation(SequenceLength, TempChromosome_1);
        end
        if (rand < MutationProbability)
            %%% Mutation by assigning a new random value
            TempChromosome_2 = random_value_mutation(SequenceLength, TempChromosome_2);
            
            %%% Mutation by order changing (also uncomment line 122 and comment out lines 119 and 126) 
            %%%     TempChromosome_2 = order_changing_mutation(SequenceLength, TempChromosome_2);
        end
        
        %%% Add to the new population if there is capacity
        PopulationNewSize = PopulationNewSize + 1;
        PopulationNew(PopulationNewSize,:) = TempChromosome_1;
        
        if (PopulationNewSize < PopulationSize)
            PopulationNewSize = PopulationNewSize + 1;
            PopulationNew(PopulationNewSize, :) = TempChromosome_2;
        end
    end
    
    Population(:, 1:SequenceLength) = PopulationNew;
end

%%% Evaluate final fitness scores and rank them
for i = 1:PopulationSize
    [Fitness, Path] = simulate_ant(Environment, Population(i, 1:SequenceLength));
    Population(i, SequenceLength+1) = Fitness;
end

Population = sortrows(Population, SequenceLength+1);

%%% Display the best (highest) fitness value
[MaxFitness, BestPath] = simulate_ant(Environment, Population(PopulationSize, 1:SequenceLength));
fprintf('Best solution: Fitness = %f\n', MaxFitness);

Map = load('muir_world.txt');
figure(3);
imagesc(Map);

for i = 1:length(BestPath)
    Map(BestPath(i, 1), BestPath(i, 2)) = -1;
end

figure(4);
imagesc(Map);

%%% Mutation by randomly selecting two mutation points and assigning new
%%% values based on which digit in the instruction set the mutation point
%%% is (can be 1-4 or 0-9)
function [TempChromosome] = random_value_mutation(SequenceLength, TempChromosome)
    MutationPoint_1 = randi([1 SequenceLength]);
    MutationPoint_1_Rem = rem(MutationPoint_1, 3);

    if MutationPoint_1_Rem ~= 1
        TempChromosome(MutationPoint_1) = randi([0 9]);
    else
        TempChromosome(MutationPoint_1) = randi([1 4]);
    end

    MutationPoint_2 = randi([1 SequenceLength]);
    MutationPoint_2_Rem = rem(MutationPoint_2, 3);

    if MutationPoint_2_Rem ~= 1
        TempChromosome(MutationPoint_2) = randi([0 9]);
    else
        TempChromosome(MutationPoint_2) = randi([1 4]);
    end
end

%%% Mutation by randomly selecting one mutation point and then randomly
%%% selecting another based on which digit in the instruction set the 
%%% mutation point is and exchanging the two values
function [TempChromosome] = order_changing_mutation(SequenceLength, TempChromosome)
    MutationPoint_1 = randi([1 SequenceLength]);
    MutationPoint_1_Rem = rem(MutationPoint_1, 3);
    
    MutationValue_1 = TempChromosome(MutationPoint_1);

    %%% Ensures that digit is only swapped with corresponding digit from another instruction set
    if MutationPoint_1_Rem == 2
        MutationPoint_2 = (randi([1 (SequenceLength/3)])*3)-2;

        TempChromosome(MutationPoint_1) = TempChromosome(MutationPoint_2);
        TempChromosome(MutationPoint_2) = MutationValue_1;
    elseif MutationPoint_1_Rem == 1
        MutationPoint_2 = (randi([1 (SequenceLength/3)])*3)-1;

        TempChromosome(MutationPoint_1) = TempChromosome(MutationPoint_2);
        TempChromosome(MutationPoint_2) = MutationValue_1;
    elseif MutationPoint_1_Rem == 0
        MutationPoint_2 = (randi([1 (SequenceLength/3)])*3);

        TempChromosome(MutationPoint_1) = TempChromosome(MutationPoint_2);
        TempChromosome(MutationPoint_2) = MutationValue_1;
    end
end
