/*
 * (1+1)-ES.h
 *
 *  Created on: 03 apr 2017
 *      Author: RedShy
 */

#ifndef SRC__1_1__ES_H_
#define SRC__1_1__ES_H_

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

int evolutionStrategy_one_one_rs(const std::vector<unsigned>& s1,
		const std::vector<unsigned>& s2, const size_t& s1l, const size_t& s2l,

		const std::vector<unsigned>& sig1, const std::vector<unsigned>& sig2,
		const size_t& sig1l, const size_t& sig2l,

		const size_t& p1, matching_schema<bool>& m, edit_distance& e,

		const unsigned max_generations)
{
	unsigned generation = 0;
	unsigned plateu = 0;
//	unsigned metaPlateu = 0;
	unsigned lastBest = std::numeric_limits<unsigned int>::max();
	const unsigned threshold = 5;
	const unsigned lapCheckPoint = 500;
	const unsigned maxPlateu = 10;



	ES_MatchingSchema parent(sig1, sig2);

	//Random start
	parent.shuffle();
	parent.costValue = e.edit_distance_matching_schema_enhanced(s1, s2, s1l,
			s2l, parent.sigma1, parent.sigma2, sig1l, sig2l, m);

	ES_MatchingSchema best = parent;
	while (generation <= max_generations)
	{
		//Produce child
		ES_MatchingSchema child = parent;

		//mutate child
		child.mutate();

//		unsigned d = e.edit_distance_matching_schema_enhanced(s1, s2, s1l, s2l,
//				child.sigma1, child.sigma2, sig1l, sig2l, m);
//		std::cout << " REALchildValue= " << d << " parentValue= "
//				<< parent.costValue << " bestValue= " << best.costValue << endl;

		//validate child
		if (ES_isValid(child))
		{
//			std::cout << " parentValueBEFOREEDIT= " << parent.costValue;
			int newDistance =
					e.edit_distance_matching_schema_enhanced_with_diagonal(s1,
							s2, s1l, s2l, child.sigma1, child.sigma2, sig1l,
							sig2l, m, parent.costValue);

//			std::cout << " newDistance= " << newDistance << endl;
			if (newDistance != -1)
			{
				//The child is better than its father, so he become new parent
				parent = child;
				parent.costValue = newDistance;

//				std::cout << "AFTER ASSIGNENT" << " childValue= "
//						<< child.costValue << " parentValue= "
//						<< parent.costValue << " bestValue= " << best.costValue
//						<< endl << endl;

				plateu = 0;

				//TODO maybe we can do better
				if (parent.costValue < best.costValue)
				{
					best = parent;
//					metaPlateu = 0;
				}
			}
			else
			{
				plateu++;
//				std::cout << "plateu= " << plateu << " WORSE!"
//						<< " childValue= " << child.costValue
//						<< " parentValue= " << parent.costValue
//						<< " bestValue= " << best.costValue << endl;
				if (plateu == maxPlateu)
				{
					plateu = 0;

//					std::cout << "BEFORE SHUFFLE " << "costValue="
//							<< parent.costValue << endl;
//					parent.print();
					parent.shuffle();
					parent.costValue = e.edit_distance_matching_schema_enhanced(
							s1, s2, s1l, s2l, parent.sigma1, parent.sigma2,
							sig1l, sig2l, m);
					if (parent.costValue < best.costValue)
					{
						best = parent;
//						metaPlateu = 0;
					}

//					std::cout << "AFTER SHUFFLE" << "costValue="
//							<< parent.costValue << endl;
//					parent.print();

				}
			}
			//else the child is worse than its father so he is discarded
		}
		else
		{
			//TODO: not valid, maybe mutate until is valid?
			//repeat iteration
			continue;
		}

		generation++;

		//TODO experimental return on investiment
		if (generation % lapCheckPoint == 0)
		{
			//few rendiment in the last checkpoint, so stop here
			if (lastBest - best.costValue <= threshold)
			{
				break;
			}

			//new value for checkpoint
			lastBest = best.costValue;
		}

//		metaPlateu++;
//		if (metaPlateu == 1000)
//		{
//			break;
//		}
	}

	//TODO return best of all
	return best.costValue;
}




#endif /* SRC__1_1__ES_H_ */
