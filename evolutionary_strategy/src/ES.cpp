#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <queue>
using namespace std;

struct MatchingSchema
{
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
			unsigned i_1 = rand() % sigma_1;
			unsigned i_2 = rand() % sigma_1;

			unsigned temp = sigma_1[i_1];
			sigma_1[i_1] = sigma_1[i_2];
			sigma_1[i_2] = temp;

			//second vector
			i_1 = rand() % sigma_2;
			i_2 = rand() % sigma_2;

			temp = sigma_2[i_1];
			sigma_2[i_1] = sigma_2[i_2];
			sigma_2[i_2] = temp;
		}

		void calculateCost()
		{
			//TODO calculate edit distance for this matching schema, maybe has to be done outside of this class?

			costValue = rand() % 999999999999;
		}

		bool operator<(const MatchingSchema& m) const
		{
			return this->costValue >= m.costValue;
		}

};

bool isValid(MatchingSchema m);
vector<MatchingSchema> selectBestIndividuals(unsigned mu,
		const vector<MatchingSchema>& parents,
		const vector<MatchingSchema>& children);
vector<MatchingSchema> selectBestIndividuals(unsigned mu,
		const vector<MatchingSchema>& children);

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
			push_heap(parents.begin(), parents.end());
		}
		else
		{
			//TODO: not valid, maybe mutate until is valid?
			//repeat iteration
			i--;
		}
	}

	vector<MatchingSchema> children;
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

				children.push_back(child);
				push_heap(children.begin(), children.end());
			}
			else
			{
				//TODO: not valid, maybe mutate until is valid?
				//repeat iteration
				i--;
			}
		}

		if (plusSelection)
		{
			//Select mu parents for next generation from parents+children
			parents = selectBestIndividuals(mu, parents, children);
		}
		else
		{
			//Comma selection
			//Select mu parents for next generation from children
			parents = selectBestIndividuals(mu, children);
		}

		children.clear();
		generation++;
	}

	//TODO return best of all
	return parents.front().costValue;
}

bool isValid(MatchingSchema m)
{
	//TODO validate a matching schema
	return true;
}

vector<MatchingSchema> selectBestIndividuals(unsigned mu,
		vector<MatchingSchema>& parents, vector<MatchingSchema>& children)
{
	//TODO using a heap

	//merge the two heaps, silly algorithm
	while (parents.size() != 0)
	{
		children.push_back(parents.back());
		parents.pop_back();
	}

	make_heap(children.begin(), children.end());

	return selectBestIndividuals(mu, children);

}
vector<MatchingSchema> selectBestIndividuals(unsigned mu,
		vector<MatchingSchema>& children)
{
	//TODO using a heap

	vector<MatchingSchema> bestIndividuals;

	for (unsigned i = 0; i < mu; i++)
	{
		pop_heap(children.begin(), children.end());
		bestIndividuals.push_back(children.back());
		children.pop_back();
	}

	return bestIndividuals;
}

int main()
{
	srand(unsigned(time(0)));
	return 0;
}
