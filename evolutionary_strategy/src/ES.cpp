#include<iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
using namespace std;

struct MatchingSchema
{
		vector<unsigned> sigma_1;
		vector<unsigned> sigma_2;
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


};

bool isValid(MatchingSchema m);
vector<MatchingSchema> selectBestIndividuals(unsigned mu,
		const vector<MatchingSchema>& parents,
		const vector<MatchingSchema>& children);
vector<MatchingSchema> selectBestIndividuals(unsigned mu,
		const vector<MatchingSchema>& parents);

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
		parents.push_back(startingMatchingSchema);
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
				children.push_back(child);
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
			//Select from parents+children
			parents = selectBestIndividuals(mu, parents, children);
		}
		else
		{
			//Comma selection
			//Select from children
			parents = selectBestIndividuals(mu, children);
		}

		children.clear();
		generation++;
	}

	//TODO return best of all
	return 0;
}

bool isValid(MatchingSchema m)
{
	//TODO validate a matching schema
	return true;
}

vector<MatchingSchema> selectBestIndividuals(unsigned mu,
		const vector<MatchingSchema>& parents,
		const vector<MatchingSchema>& children)
{
	//TODO
}
vector<MatchingSchema> selectBestIndividuals(unsigned mu,
		const vector<MatchingSchema>& parents)
{
	//TODO
}


int main()
{
	srand(unsigned(time(0)));
	return 0;
}
