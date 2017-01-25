#!/usr/bin/python

# dumb idea: selecting pairs two by two - correct, but super slow
from random import random, sample, choice


N = 100 # population size
n = 10  # number of people with bad allele
N_ITER = 100 # number of simulation iterations

number_of_generations = 100

class Individual:
    def __init__ (self, genotype):
        self.parent    = None
        self.children  = []
        self.genotype  = genotype

###############
def initial_population():
    population = []
    for i in range(N):
        # x carry bad gene
        population.append(Individual(["G","G"]))

    for carrier in sample(population, n):
        carrier.genotype = ["G","B"]
    return population

##########
def procreate(current_generation):
    next_generation = []

    while len(current_generation) > 1:
        # they mate at random
        parents = sample(current_generation, 2)
        for parent in parents:
            current_generation.remove(parent)
        for i in range(2):
            genotype = [choice(parents[0].genotype), choice(parents[1].genotype)]
            next_generation.append(Individual(genotype))
    return next_generation

####################################
def find_allele_fraction(population):
    bcount = 0.0
    for x in population:
        for a in x.genotype:
            if a=="B": bcount +=1

    return bcount/(len(population)*2) # two alleles per individual

####################################
def main():


    avg_homozygous_ct   = [0]*number_of_generations
    avg_allele_fraction = [0]*number_of_generations

    for it in range(N_ITER):

        generation = [[],[]]
        current_gen = 0
        next_gen = 1
        generation[0] = initial_population()

        for gen_number in range(number_of_generations):
            # how many are homozygous as a function of time (number of generations)
            homozygous_ct = len([x for x in generation[current_gen] if set(x.genotype)=={"B","B"}])
            allele_fraction = find_allele_fraction(generation[current_gen])
            #print " %3d  %5.3f " %(gen_number, allele_fraction)
            avg_homozygous_ct[gen_number] += homozygous_ct
            avg_allele_fraction[gen_number] += allele_fraction
            generation[next_gen] = procreate(generation[current_gen])
            current_gen = 1 - current_gen
            next_gen = 1 - current_gen
        #exit(1)
    for gen_number in range(number_of_generations):
        avg_homozygous_ct[gen_number] /= N_ITER*1.0
        avg_allele_fraction[gen_number] /= N_ITER*1.0

    for gen_number in range(number_of_generations):
        print " %3d  %5.1f   %5.3f " %(gen_number, avg_homozygous_ct[gen_number], avg_allele_fraction[gen_number])


####################################
if __name__ == '__main__':
    main()
