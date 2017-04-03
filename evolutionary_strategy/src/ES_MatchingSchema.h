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
				costValue(0), sigma1l(
						_sigma_1.size()), sigma2l(_sigma_2.size())
		{
			sigma1 = new unsigned[_sigma_1.size()];
			std::iota(sigma1, sigma1 + _sigma_1.size(), 0);

			sigma2 = new unsigned[_sigma_2.size()];
			std::iota(sigma2, sigma2 + _sigma_2.size(), 0);
		}

		ES_MatchingSchema(const ES_MatchingSchema& m) :
				sigma1l(m.sigma1l), sigma2l(m.sigma2l), costValue(m.costValue)
		{
			sigma1 = new unsigned[sigma1l];
			std::copy(m.sigma1, m.sigma1 + sigma1l, sigma1);

			sigma2 = new unsigned[sigma2l];
			std::copy(m.sigma2, m.sigma2 + sigma2l, sigma2);
		}

		void shuffle()
		{
			std::random_shuffle(sigma1, sigma1 + sigma1l);
			std::random_shuffle(sigma2, sigma2 + sigma2l);
		}

		void mutate()
		{
			//Perform a single, simple swap of two indices for every vector

			//first vector
			unsigned i_1 = rand() % sigma1l;
			unsigned i_2 = rand() % sigma1l;

			unsigned temp = sigma1[i_1];
			sigma1[i_1] = sigma1[i_2];
			sigma1[i_2] = temp;

			//second vector
			i_1 = rand() % sigma2l;
			i_2 = rand() % sigma2l;

			temp = sigma2[i_1];
			sigma2[i_1] = sigma2[i_2];
			sigma2[i_2] = temp;
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

		ES_MatchingSchema& operator=(const ES_MatchingSchema& m)
		{
			costValue = m.costValue;

			//TODO: we can assume that every child has the same simgas length
//			if(sigma1l!=m.sigma1l)
//			{
//				sigma1l = m.sigma1l;
//				delete[] sigma1;
//				sigma1 = new unsigned[sigma1l];
//			}
			std::copy(m.sigma1, m.sigma1 + sigma1l, sigma1);

//			if(sigma2l!=m.sigma2l)
//			{
//				sigma2l = m.sigma2l;
//				delete[] sigma2;
//				sigma2 = new unsigned[sigma2l];
//			}
			std::copy(m.sigma2, m.sigma2 + sigma2l, sigma2);
			return *this;
		}

		~ES_MatchingSchema()
		{
			delete[] sigma1;
			delete[] sigma2;
		}



		unsigned* sigma1;
		const size_t sigma1l;

		unsigned* sigma2;
		const size_t sigma2l;

		unsigned costValue;
};



#endif /* SRC_ES_MATCHINGSCHEMA_H_ */
