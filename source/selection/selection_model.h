#pragma once

#include <stdio.h>
#include <vector>
#include <map>
#include <string>


using namespace std;


namespace mimir {
    
	struct Sequence { 
    public:
			string aminoacide;
			string V;
			string J;
			int* aminoacide_indexes;
			int Vindex;
			int Jindex;
			bool valid;
			Sequence(string aa,string v,string j):aminoacide(aa),V(v),J(j)
			{
				valid=true;
				aminoacide_indexes=NULL;
			};
        
			~Sequence(){
				delete[] aminoacide_indexes;
				
			}
			void convertToIndexes(map<string,int>& V_indexes, map<string,int>& J_indexes, map<char,int>& AA_indexes){
				aminoacide_indexes=new int[aminoacide.length()];
				for(int i=0;i<aminoacide.length();i++){
					if(aminoacide[i]!='*'&&aminoacide[i]!='~')
						aminoacide_indexes[i]=(AA_indexes)[aminoacide[i]];
					else{
						valid=false;
						break;
					}
				}
				Vindex=V_indexes[V];
				Jindex=J_indexes[J];
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
		void fit(const SequenceVector &data_seq,  SequenceVector &gen_seq);


        /**
         *
         */
		float predict(const Sequence &seq);
		float* predictMany(const SequenceVector &seq);


		float* get_q_L(){return q_L;};
		float* get_q_VJ(){
			return q_VJ;
		};
		float* get_ilA(){return q_ilA;};
		float getZ(){return Z;};
		int getMinL(){return minL;};
		int getMaxL(){return maxL;};

		map<string,int>* getVIndexes(){return V_indexes;};
		map<string,int>* getJIndexes(){return J_indexes;};
		map<char,int>* getAAindexes(){return &aminoAcidIndexes;};

        
	private:
		static const float EPS;
		static const int MAX_STEP;
		int minL,maxL;
		map<string,int>* V_indexes, * J_indexes;
		map<char,int> SelectionModel::aminoAcidIndexes;

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

		void evalf_gen_Ldistribution(const SequenceVector &gen_seq, float* l_distribution);
		void evalf_gen_VJdistribution(const SequenceVector &gen_seq,float* VJ_distribution);
		void evalf_gen_AAdistribution(const SequenceVector &gen_seq,float* AA_distribution);
		float evalf_Z(const SequenceVector &gen_seq);

		void transformData(SequenceVector *seq);
	};
}

