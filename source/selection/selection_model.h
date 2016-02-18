#pragma once

#include <stdio.h>
#include <vector>
#include <map>
#include <string>


using namespace std;


namespace mimir {
    
	struct Sequence { 
    public:
			string* aminoacide;
			string* V;
			string* J;
			Sequence(string* aa,string* v,string* j):aminoacide(aa),V(v),J(j)
			{};
        
			~Sequence(){
				delete aminoacide;
				delete V;
				delete J;
			}
        
    };


    typedef std::vector<Sequence> SequenceVector;


	class SelectionModel {
	public:

        /**
         * DOCS HERE
         */
		SelectionModel(void);


        /**
         *
         */
		~SelectionModel(void);


        /**
         *
         */
		void fit(const SequenceVector &data_seq, const SequenceVector &gen_seq);


        /**
         *
         */
		void predict(const Sequence &seq);
        
	private:
		static const float EPS;
		static const int MAX_STEP;
		int minL,maxL;
		map<string,int>* V_indexes, * J_indexes;


		float* data_Ldistribution;
		float* data_VJpairDistribution;
		float* data_AAdistibution;

		float* q_L;
		float* q_VJ;
		float* q_ilA;
		float Z;

		float Q(const Sequence &seq);

		inline float getLProbabilityInData(int L){
			if(L<minL||L>maxL)
				return 0;
			return data_Ldistribution[L];
		}
		void findMinMaxLength(const SequenceVector &data_seq, const SequenceVector &gen_seq);
		float* evalfDataLDistribution(const SequenceVector &data_seq, int minL, int maxL, int minFrequency);
		map<string,int>* extractVSet(const SequenceVector &data_seq, const SequenceVector &gen_seq);
		map<string,int>* extractJSet(const SequenceVector &data_seq, const SequenceVector &gen_seq);

		void evalf_gen_Ldistribution(const SequenceVector &gen_seq, int minFrequency, float* l_distribution);
		float evalf_Z(const SequenceVector &gen_seq);

	};
}

