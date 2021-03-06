#!/usr/bin/python
from random import shuffle, sample, choice

N = 100 # population size
n = 10  # number of people with bad allele
N_ITER = 100 # number of simulation iterations


number_of_generations = 100
nr_of_generations_to_track = 3

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
def procreate(generation):
    next_generation = []
    current_generation = generation[0][:] # this is a list copy
    shuffle(current_generation)
    for p in xrange(0,len(current_generation),2):
        parent1 = current_generation[p]
        parent2 = current_generation[p+1]
        for i in range(2):
            genotype = [choice(parent1.genotype), choice(parent2.genotype)]
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

        generation  = [[]]*nr_of_generations_to_track
        current_gen = 0
        generation[0] = initial_population()

        for gen_number in range(number_of_generations):
            # how many are homozygous as a function of time (number of generations)
            homozygous_ct = len([x for x in generation[current_gen] if set(x.genotype)=={"B","B"}])
            allele_fraction = find_allele_fraction(generation[current_gen])
            #print " %3d  %5.3f " %(gen_number, allele_fraction)
            avg_homozygous_ct[gen_number] += homozygous_ct
            avg_allele_fraction[gen_number] += allele_fraction
            next_generation = procreate(generation)
            for g in reversed(range(nr_of_generations_to_track-1)):
                generation[g+1] = generation[g]
            generation[0] = next_generation
        #exit(1)

    for gen_number in range(number_of_generations):
        avg_homozygous_ct[gen_number] /= N_ITER*1.0
        avg_allele_fraction[gen_number] /= N_ITER*1.0

    for gen_number in range(number_of_generations):
        print " %3d  %5.1f   %5.3f " %(gen_number, avg_homozygous_ct[gen_number], avg_allele_fraction[gen_number])


####################################
if __name__ == '__main__':
    main()
