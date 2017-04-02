#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <queue>
using namespace std;

unsigned n = 0;
struct MatchingSchema
{
		MatchingSchema(const vector<unsigned>&_sigma_1,
				const vector<unsigned>&_sigma_2) :
				sigma_1(_sigma_1), sigma_2(_sigma_2)
		{

		}

		vector<unsigned> sigma_1;
		vector<unsigned> sigma_2;
		unsigned costValue;
		void shuffle()
		{
			random_shuffle(sigma_1.begin(), sigma_1.end());
			random_shuffle(sigma_2.begin(), sigma_2.end());
		}

		void mutate()
		{
			//Perform a single, simple swap of two indices for every vector

			//first vector
			unsigned i_1 = rand() % sigma_1.size();
			unsigned i_2 = rand() % sigma_1.size();

			unsigned temp = sigma_1[i_1];
			sigma_1[i_1] = sigma_1[i_2];
			sigma_1[i_2] = temp;

			//second vector
			i_1 = rand() % sigma_2.size();
			i_2 = rand() % sigma_2.size();

			temp = sigma_2[i_1];
			sigma_2[i_1] = sigma_2[i_2];
			sigma_2[i_2] = temp;
		}

		void calculateCost()
		{
			//TODO calculate edit distance for this matching schema, maybe has to be done outside of this class?

			costValue = rand() % 999999999;
			n++;
		}

		bool operator<(const MatchingSchema& m) const
		{
			return this->costValue >= m.costValue;
		}

};

bool isValid(MatchingSchema m);
vector<MatchingSchema> selectBestIndividuals(const unsigned mu,
		vector<MatchingSchema>& individuals, const bool plusSelection);

unsigned evolutionStrategy(const unsigned max_generations, const unsigned mu,
		const unsigned lambda, const bool plusSelection,
		MatchingSchema startingMatchingSchema)
{
	unsigned generation = 0;

	//Generate mu random individuals
	vector<MatchingSchema> parents;
	for (unsigned i = 0; i < mu; i++)
	{
		startingMatchingSchema.shuffle();

		//validate matching schema
		if (isValid(startingMatchingSchema))
		{
			startingMatchingSchema.calculateCost();

			parents.push_back(startingMatchingSchema);
//			push_heap(parents.begin(), parents.end());
		}
		else
		{
			//TODO: not valid, maybe mutate until is valid?
			//repeat iteration
			i--;
		}
	}

	while (generation <= max_generations)
	{
		//Generate lambda children. Only mutation, no recombination
		for (unsigned i = 0; i < lambda; i++)
		{
			//Choose random parent
			unsigned p = rand() % mu;

			//Produce child, in the case parents=1 (like this) just clone
			MatchingSchema child = parents[p];

			//mutate child
			child.mutate();

			//validate child
			if (isValid(child))
			{
				child.calculateCost();

				parents.push_back(child);
//				push_heap(children.begin(), children.end());
			}
			else
			{
				//TODO: not valid, maybe mutate until is valid?
				//repeat iteration
				i--;
			}
		}


		parents = selectBestIndividuals(mu, parents, plusSelection);


		generation++;
	}

	//TODO return best of all
	make_heap(parents.begin(), parents.end());
	return parents.front().costValue;
}

bool isValid(MatchingSchema m)
{
	//TODO validate a matching schema
	return true;
}

vector<MatchingSchema> selectBestIndividuals(const unsigned mu,
		vector<MatchingSchema>& individuals, const bool plusSelection)
{
	//TODO using a heap
	vector<MatchingSchema> bestIndividuals;
	if (plusSelection)
	{
		make_heap(individuals.begin(), individuals.end());

		for (unsigned i = 0; i < mu; i++)
		{
			pop_heap(individuals.begin(), individuals.end());
			bestIndividuals.push_back(individuals.back());
			individuals.pop_back();
		}
	}
	else
	{
		//start from mu instead of 0
		vector<MatchingSchema>::iterator it = individuals.begin();
		it += mu;
		make_heap(it, individuals.end());

		for (unsigned i = 0; i < mu; i++)
		{
			pop_heap(it, individuals.end());
			bestIndividuals.push_back(individuals.back());
			individuals.pop_back();
		}
	}

	return bestIndividuals;
}

int main()
{
	srand(unsigned(time(0)));
	unsigned A1 = 500;
	unsigned A2 = 500;
	vector<unsigned> s1;
	for (unsigned i = 0; i < A1; i++)
	{
		s1.push_back(i);
	}
	vector<unsigned> s2;
	for (unsigned i = 0; i < A2; i++)
	{
		s2.push_back(i);
	}
	MatchingSchema m1(s1, s2);
	clock_t begin = clock();

	cout << evolutionStrategy(10000, 50, 50, true, m1) << endl;

	clock_t end = clock();
	cout << double(end - begin) / CLOCKS_PER_SEC << endl;
	cout << n << endl;
	return 0;
}
