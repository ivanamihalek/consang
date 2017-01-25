#!/usr/bin/python
from random import shuffle, sample, choice
from sys import argv

N = 500 # population size - must be even
n = 10 # number of people with bad allele
N_ITER = 300 # number of simulation iterations

number_of_generations = 100
nr_of_generations_to_track = 1

###############
class Individual:
    def __init__ (self, genotype):
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
def set_children(current_generation, next_generation, index1, index2):
    parent1 = current_generation[index1]
    parent2 = current_generation[index2]
    children = []
    for i in range(2):
        genotype = [choice(parent1.genotype), choice(parent2.genotype)]
        child = Individual(genotype)
        children.append(child)
    next_generation[index1] = children[0]
    next_generation[index2] = children[1]


##########
def random_mate(current_generation):
    next_generation = [None]*N
    shuffle(current_generation)
    for p in xrange(0,len(current_generation),2):
        set_children(current_generation, next_generation,p, p+1)
    return next_generation


##########
def consanguinous_mate(current_generation, offset):
    next_generation = [None]*N
    for p in xrange(offset,len(current_generation),2):
         set_children(current_generation, next_generation, p, (p+1)%N)

    return next_generation

##########
def procreate(current_generation, gen_number):
    next_generation = random_mate(current_generation)
    #next_generation = consanguinous_mate(current_generation, gen_number%2)
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
        if not it%10: print "iteration:", it
        generation = [[]]*nr_of_generations_to_track
        # need three generations to figure out cousins
        generation[0] = initial_population()
        current_gen = 0

        for gen_number in range(number_of_generations):
            # how many are homozygous as a function of time (number of generations)
            homozygous_ct = len([x for x in generation[current_gen] if x.genotype==["B","B"]])
            allele_fraction = find_allele_fraction(generation[current_gen])
            #cf = consang_parent_fraction(generation[current_gen])
            #print " gen number: %3d   allele fraction: %5.3f  homozygous: %3d   cons parent fraction:  %5.3f" \
            #      %(gen_number, allele_fraction, homozygous_ct, cf)
            avg_homozygous_ct[gen_number] += homozygous_ct
            avg_allele_fraction[gen_number] += allele_fraction
            next_generation = procreate(generation[0], gen_number)
            generation[0] = next_generation

    for gen_number in range(number_of_generations):
        avg_homozygous_ct[gen_number] /= N_ITER*1.0
        avg_allele_fraction[gen_number] /= N_ITER*1.0

    for gen_number in range(number_of_generations):
        print " %3d  %5.1f   %5.3f " %(gen_number, avg_homozygous_ct[gen_number], avg_allele_fraction[gen_number])

    print "N:",N
    print "n:",n

####################################
if __name__ == '__main__':
    main()
