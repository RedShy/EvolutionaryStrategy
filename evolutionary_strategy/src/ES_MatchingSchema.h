/*
 * ES_MatchingSchema.h
 *
 *  Created on: 02 apr 2017
 *      Author: RedShy
 */
#ifndef SRC_ES_MATCHINGSCHEMA_H_
#define SRC_ES_MATCHINGSCHEMA_H_
#include<vector>

struct ES_MatchingSchema
{
		ES_MatchingSchema(const std::vector<unsigned>&_sigma_1,
				const std::vector<unsigned>&_sigma_2) :
				sigma_1(_sigma_1), sigma_2(_sigma_2), costValue(0)
		{

		}

		std::vector<unsigned> sigma_1;
		std::vector<unsigned> sigma_2;
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
		}

		bool operator<(const ES_MatchingSchema& m) const
		{
			return this->costValue >= m.costValue;
		}

};



#endif /* SRC_ES_MATCHINGSCHEMA_H_ */
