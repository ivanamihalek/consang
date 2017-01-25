#!/usr/bin/python

# dumb for  two reasons
# 1) sampling from generation is small
# 2) cousins should come from families that have tradition of marrying cousins
#    not from  a different family in each generation

from random import random, sample, choice
from sys import argv

N = 500 # population size
n = 10 # number of people with bad allele
N_ITER = 100 # number of simulation iterations
consang_fraction = 0.5

number_of_generations = 100
nr_of_generations_to_track = 3

###############
class Individual:
    def __init__ (self, genotype):
        self.parents   = []
        self.children  = []
        self.married   = False
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
def random_mate(current_generation):
    parents = sample(current_generation, 2)
    return parents

##########
def find_cousins(individual):
    grandparents = set()
    siblings = set()
    for parent in individual.parents:
        grandparents |= set(parent.parents)
        siblings |= set(parent.children)
    grampas_grandchildren = set()
    for grampa in grandparents:
        for uncle in grampa.children:
            grampas_grandchildren |= set(uncle.children)
    cousins = grampas_grandchildren - siblings
    return cousins

##########
def consanguinous_mate(current_generation):
    parents = [] # we will return empty list if there are no more cousins available
    while len(current_generation) > 1 and parents==[]:
        parent1 = choice(current_generation)
        cousins = find_cousins(parent1)
        unmarried_cousins = [x for x in cousins if not x.married]
        if len(unmarried_cousins) == 0: continue
        parent2 = choice(unmarried_cousins)
        parents = [parent1, parent2]
    #if no cousins available, we return an empty parents list
    return parents

##########
def procreate(population, consang_fraction = 0):
    next_generation = []*nr_of_generations_to_track
    current_generation = population[:] # this is a list copy
    fraction_married  = 0.0
    done_with_cousins = fraction_married >= consang_fraction
    while len(current_generation) > 1:
        # marry the cousins first, so we do not run out of choices
        if not done_with_cousins:
            parents = consanguinous_mate(current_generation)
            fraction_married = float(N-len(current_generation))/N
            done_with_cousins = (parents ==[]) or fraction_married >= consang_fraction

        else:
            parents = random_mate(current_generation)

        for parent in parents:
            parent.married = True
            current_generation.remove(parent)

        for i in range(2):
            genotype = [choice(parents[0].genotype), choice(parents[1].genotype)]
            child = Individual(genotype)
            [parent.children.append(child) for parent in parents]
            child.parents = parents
            next_generation.append(child)

    return next_generation

####################################
def find_allele_fraction(population):
    bcount = 0.0
    for x in population:
        for a in x.genotype:
            if a=="B": bcount +=1

    return bcount/(len(population)*2) # two alleles per individual
####################################
def consang_parent_fraction(population):
    cp_count = 0.0
    for x in population:
        cousins = find_cousins(x.parents[0])
        if x.parents[1] in cousins: cp_count += 1
    return cp_count/len(population)


####################################
def main():
    global consang_fraction # needed to modify global variable
    if len(argv)>1: consang_fraction = float(argv[1])
    avg_homozygous_ct   = [0]*number_of_generations
    avg_allele_fraction = [0]*number_of_generations

    for it in range(N_ITER):
        if not it%10: print "iteration:", it
        generation = [[]]*nr_of_generations_to_track
        # need three generations to figure out cousins
        generation[2] = initial_population()
        generation[1] = procreate(generation[2])
        generation[0] = procreate(generation[1])
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
            next_generation = procreate(generation[0], consang_fraction = consang_fraction)
            for g in reversed(range(nr_of_generations_to_track-1)):
                generation[g+1] = generation[g]
            generation[0] = next_generation

    for gen_number in range(number_of_generations):
        avg_homozygous_ct[gen_number] /= N_ITER*1.0
        avg_allele_fraction[gen_number] /= N_ITER*1.0

    for gen_number in range(number_of_generations):
        print " %3d  %5.1f   %5.3f " %(gen_number, avg_homozygous_ct[gen_number], avg_allele_fraction[gen_number])

    print "N:",N
    print "n:",n
    print "consang fraction:", consang_fraction

####################################
if __name__ == '__main__':
    main()
