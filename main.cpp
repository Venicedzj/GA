#include<iostream>
#include<vector>
#include <random>
using namespace std;

#define POPULATION_NUM 100
#define DIMENSION 100
#define Pc 0.7
#define Pm 0.07
#define UB 5.12
#define LB -5.12
#define GEN 3000
int gen = 1;

struct Individual {
	vector<double> gene;
	double fitness;
};

void Init();
double CalFitness(vector<double>& gene);
void Print();
void Selection();
void Crossover();
void Mutation();
void Evaluation();
void Elitist();
inline double random(double a, double b) { return ((double)rand() / RAND_MAX) * (b - a) + a; }
inline int random_int(int a, int b) { return rand() % (b - a) + a; }

vector<Individual> population(POPULATION_NUM);
vector<Individual> offspring(POPULATION_NUM);
Individual best_ind;

int main() {
	srand((int)time(NULL));
	Init();
	Print();
	while (gen++ < GEN) {
		Selection();
		Crossover();
		Mutation();
		Evaluation();
		Elitist();
		Print();
	}
}

/* initialize population with random number between lower bound(LB) and upper bound(UB) */
void Init() {
	best_ind.fitness = DBL_MAX;
	for (int i = 0; i < POPULATION_NUM; ++i) {
		for (int j = 0; j < DIMENSION; ++j) {
			population[i].gene.push_back(random(LB, UB));
		}
		population[i].fitness = CalFitness(population[i].gene);

		//get the best individual's index
		if (population[i].fitness < best_ind.fitness) {
			best_ind = population[i];
		}
	}
}

void Selection() {
	vector<double> rate_intersc(POPULATION_NUM);
	double sum_fitness = 0.0;
	for (auto i : population) {
		sum_fitness += 1.0 / i.fitness;
	}

	rate_intersc[0] = (1.0 / population[0].fitness) / sum_fitness;
	for (int i = 1; i < POPULATION_NUM; ++i) {
		rate_intersc[i] = rate_intersc[i - 1] + ((1.0 / population[i].fitness) / sum_fitness);
	}

	for (int i = 0; i < POPULATION_NUM; ++i) {
		double rnd = random(0.0, 1.0);
		int k = 0;
		while (rnd > rate_intersc[k] && k < POPULATION_NUM - 1) { k++; }
		offspring[i] = population[k];
	}
	population = offspring;
}

void Crossover() {
	vector<int> corss_pool;
	for (int i = 0; i < POPULATION_NUM; ++i) {
		double rnd_Pc = random(0.0, 1.0);
		if (rnd_Pc < Pc) {
			corss_pool.push_back(i);
			if (corss_pool.size() == 2) {
				int cr_point = random_int(0, DIMENSION);
				for (int j = 0; j < cr_point; ++j)
					swap(population[corss_pool[0]].gene[j], population[corss_pool[1]].gene[j]);
				corss_pool.clear();
			}
		}
	}
}

void Mutation() {
	for (int i = 0; i < POPULATION_NUM; ++i) {
		for (int j = 0; j < DIMENSION; ++j) {
			double rnd = random(0.0, 1.0);
			if (rnd < Pm) {
				population[i].gene[j] = random(LB, UB);
			}
		}
	}
}

void Elitist() {
	int best, worst;
	double best_fitness = DBL_MAX;
	double worst_fitness = DBL_MIN;
	for (int i = 0; i < POPULATION_NUM; ++i) {
		if (population[i].fitness < best_fitness) {
			best = i; 
			best_fitness = population[i].fitness;
		}
		if (population[i].fitness > worst_fitness) {
			worst = i;
			worst_fitness = population[i].fitness;
		}
	}
	if (best_fitness < best_ind.fitness) best_ind = population[best];
	else population[worst] = best_ind;
}

/* select the offspring with smaller fitness as new population */
void Evaluation() {
	for (int i = 0; i < POPULATION_NUM; ++i) {
		population[i].fitness = CalFitness(population[i].gene);
	}
}

/* fitness evaluation function */
double CalFitness(vector<double>& gene) {
	double fitness = 0.0;

	///Sphere function
	for (auto i : gene) {
		fitness += i * i;
	}

	///Ackley function
	//double sum1 = 0.0, sum2 = 0.0;
	//for (auto i : pos) {
	//	sum1 += pow(i, 2);
	//	sum2 += cos(2 * PI * i);
	//}
	//fitness = 20 + e - 20 * exp(-0.2 * sqrt(sum1 / D)) - exp(sum2 / D);

	return fitness;
}

/* print the individual with smallest fitness */
void Print() {
	cout << "Generation: " << gen << " best fitness: " << best_ind.fitness << endl;
}