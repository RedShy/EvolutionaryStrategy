#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <queue>
#include "ES_MatchingSchema.h"
#include "EditDistance.h"
#include "MatchingSchema.h"


bool ES_isValid(ES_MatchingSchema m);
std::vector<ES_MatchingSchema> selectBestIndividuals(const unsigned mu,
		std::vector<ES_MatchingSchema>& individuals, const bool plusSelection);

int evolutionStrategy(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations, const unsigned mu,
		const unsigned lambda, const bool plusSelection)
{
	unsigned generation = 0;

	//TODO temporay for edit distance with diagonal
	unsigned threshold = std::numeric_limits<unsigned int>::max();

	ES_MatchingSchema startingMS(sig1, sig2);

	//Generate mu random individuals
	std::vector<ES_MatchingSchema> parents;
	for (unsigned i = 0; i < mu; ++i)
	{
		startingMS.shuffle();

		//validate matching schema
		if (ES_isValid(startingMS))
		{
//			startingMS.calculateCost();

			startingMS.costValue =
					e.edit_distance_matching_schema_enhanced_with_diagonal(s1,
							s2, s1l, s2l, startingMS.sigma1, startingMS.sigma2,
							sig1l, sig2l, m, threshold);
			parents.push_back(startingMS);
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
			ES_MatchingSchema child = parents[p];

			//mutate child
			child.mutate();

			//validate child
			if (ES_isValid(child))
			{
//				child.calculateCost();

				child.costValue =
						e.edit_distance_matching_schema_enhanced_with_diagonal(
								s1, s2, s1l, s2l, child.sigma1, child.sigma2,
								sig1l, sig2l, m, threshold);

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

bool ES_isValid(ES_MatchingSchema m)
{
	//TODO validate a matching schema
	return true;
}

std::vector<ES_MatchingSchema> selectBestIndividuals(const unsigned mu,
		std::vector<ES_MatchingSchema>& individuals, const bool plusSelection)
{
	//TODO using a heap
	std::vector<ES_MatchingSchema> bestIndividuals;
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
		//TODO Comma selection
	}

	return bestIndividuals;
}

int main()
{
	srand(unsigned(time(0)));
	unsigned A1 = 500;
	unsigned A2 = 500;
	std::vector<unsigned> s1;
	for (unsigned i = 0; i < A1; i++)
	{
		s1.push_back(i);
	}
	std::vector<unsigned> s2;
	for (unsigned i = 0; i < A2; i++)
	{
		s2.push_back(i);
	}
	ES_MatchingSchema m1(s1, s2);
	clock_t begin = clock();

//	cout << evolutionStrategy(10000, 50, 50, true, m1) << endl;

	clock_t end = clock();
	std::cout << double(end - begin) / CLOCKS_PER_SEC << std::endl;
	return 0;
}
