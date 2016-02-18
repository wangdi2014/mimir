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
		void fit(vector<Sequence>* data_seq,vector<Sequence>* gen_seq);


        /**
         *
         */
		void predict(Sequence* seq);
        
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

		float Q(Sequence* seq);

		inline float getLProbabilityInData(int L){
			if(L<minL||L>maxL)
				return 0;
			return data_Ldistribution[L];
		}
		void findMinMaxLength(vector<Sequence>* data_seq,vector<Sequence>* gen_seq);
		float* evalfDataLDistribution(vector<Sequence>* data_seq,int minL,int maxL,int minFrequency);
		map<string,int>* extractVSet(vector<Sequence>* data_seq,vector<Sequence>* gen_seq);
		map<string,int>* extractJSet(vector<Sequence>* data_seq,vector<Sequence>* gen_seq);

		void evalf_gen_Ldistribution(vector<Sequence>* gen_seq,int minFrequency,float* l_distribution);
		float evalf_Z(vector<Sequence>* gen_seq);

	};
}

