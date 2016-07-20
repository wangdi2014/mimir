#pragma once

#include <stdio.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <Eigen/Dense>
#include<Eigen/Sparse>


using namespace std;


namespace mimir {
    
	struct Sequence { 
    public:
			string aminoacide;
			vector<string> V;
			vector<string> J;
			int* aminoacide_indexes;
			int* Vindex;
			int* Jindex;
			float W_count;
			float* weights;
			bool valid;
			Sequence(string aa,vector<string> v,vector<string> j):aminoacide(aa),V(v),J(j)
			{
				valid=true;
				W_count=1;
				aminoacide_indexes=NULL;
			};
			Sequence(string aa,vector<string> v,vector<string> j,float w):aminoacide(aa),V(v),J(j),W_count(w)
			{
				valid=true;
				aminoacide_indexes=NULL;
			};
			Sequence(string aa,string v,string j):aminoacide(aa)
			{
				V.push_back(v);
				J.push_back(j);
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
				Vindex=new int[V.size()];
				Jindex=new int[J.size()];
				for(int i=0;i<V.size();i++){
					Vindex[i]=V_indexes[V[i]];
				}
				for(int i=0;i<J.size();i++){
					Jindex[i]=J_indexes[J[i]];
				}
				weights=new float[V.size()*J.size()];
				for(int i=0;i<V.size()*J.size();i++){
					weights[i]=1.0/(V.size()*J.size());
				}
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
		void fit(SequenceVector &data_seq,  SequenceVector &gen_seq,int max_iter);
		void evalfP_gen(SequenceVector &gen_seq);
		double evalfMaxDelta(SequenceVector &gen_seq,double alpha,double* F,double* emptyArray,int Lsize,int VJsize,int AAsize);
		double evalfMaxLikehood(SequenceVector &gen_seq,SequenceVector &data_seq,double alpha,double* F,int Lsize,int VJsize,int AAsize);
		double optimizeStep(SequenceVector &gen_seq,SequenceVector &data_seq,double* F,double* emptyArray,int Lsize,int VJsize,int AAsize,double* maxDelt);

        /**
         *
         */
		double predict(const Sequence &seq);
		double* predictMany(const SequenceVector &seq);

		double* Gw ;
		double* Dw ;

		double* get_q_L(){return q_L;};
		double* get_q_VJ(){return q_VJ;};
		double* get_ilA(){return q_ilA;};

		double* get_cluster_L(int i) { return cluster_L+i*(maxL - minL + 1); };
		double* get_cluster_VJ(int i) { return cluster_VJ+i*(V_indexes->size()*J_indexes->size()); };
		double* get_cluster_ilA(int i) { return cluster_ilA+i*((maxL - minL + 1)*maxL*aminoAcidIndexes.size()); };

		double* get_Pgen_L(){return gen_L;};
		double* get_Pgen_VJ(){	return gen_VJ;};
		double* get_Pgen_ilA(){return gen_ilA;};
		double getZ(){return Z;};
		int getMinL(){return minL;};
		int getMaxL(){return maxL;};

		map<string,int>* getVIndexes(){return V_indexes;};
		map<string,int>* getJIndexes(){return J_indexes;};
		map<char,int>* getAAindexes(){return &aminoAcidIndexes;};

        
	private:
		static const double EPS;
		static const int MAX_STEP;
		int minL,maxL;
		map<string,int>* V_indexes, * J_indexes;
		map<char,int> SelectionModel::aminoAcidIndexes;

		double* data_Ldistribution;
		double* data_VJpairDistribution;
		double* data_AAdistibution;

		double* q_L;
		double* q_VJ;
		double* q_ilA;

		double* cluster_L;
		double* cluster_VJ;
		double* cluster_ilA;
		
		

		double* gen_L;
		double* gen_VJ;
		double* gen_ilA;



		double Z;

		double Q(const Sequence &seq);

		inline double getLProbabilityInData(int L){
			if(L<minL||L>maxL)
				return 0;
			return data_Ldistribution[L];
		}
		void findMinMaxLength(const SequenceVector &data_seq, const SequenceVector &gen_seq);
		double* evalfDataLDistribution(const SequenceVector &data_seq, int minL, int maxL, int minFrequency);
		map<string,int>* extractVSet(const SequenceVector &data_seq, const SequenceVector &gen_seq);
		map<string,int>* extractJSet(const SequenceVector &data_seq, const SequenceVector &gen_seq);

		void evalf_gen_Ldistribution(const SequenceVector &gen_seq, double* l_distribution);
		void evalf_gen_VJdistribution(const SequenceVector &gen_seq,double* VJ_distribution);
		void evalf_gen_AAdistribution(const SequenceVector &gen_seq,double* AA_distribution);
		double evalf_Z(const SequenceVector &gen_seq);

		void transformData(SequenceVector *seq);
	};
}

