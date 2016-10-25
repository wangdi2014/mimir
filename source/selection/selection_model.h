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
    
	struct Martrix2d {
	public:
		vector<double> data;
		int W, H;
		Martrix2d(int W, int H):W(W),H(H) {
			data.resize(W*H);
			fill(0);
		}
		Martrix2d(){
			W = 0;H = 0;
		}
		void fill(double val) {
			std::fill(std::begin(data), std::end(data), 0);
		}
		double get(int i, int j) {
			if (i > W)
				printf("ERRRRRROOOOOOOOOO!");
			if (j > H)
				printf("ERRRRRROOOOOOOOOO!");
			return data[H*i + j];

		}
		void set(int i, int j, double val) {
			if (i > W)
				printf("ERRRRRROOOOOOOOOO!");
			if (j > H)
				printf("ERRRRRROOOOOOOOOO!");
			data[H*i + j] = val;
		}
		void add(int i, int j, double val) {
			if (i > W)
				printf("ERRRRRROOOOOOOOOO!");
			if (j > H)
				printf("ERRRRRROOOOOOOOOO!");
			data[H*i + j] += val;
		}
	};

	struct Martrix3d {
	public:
		vector<double> data;
		int W, H, D;
		Martrix3d(int W, int H, int D) :W(W), H(H), D(D) {
			data.resize(W*H*D);
			fill(0);
		}
		Martrix3d() :W(0), H(0), D(0) {

		}
		void fill(double val) {
			std::fill(std::begin(data), std::end(data), val);
		}
		double get(int i, int j, int k) {
			if (i > W)
				printf("ERRRRRROOOOOOOOOO!");
			if (j > H)
				printf("ERRRRRROOOOOOOOOO!");
			if (k > D)
				printf("ERRRRRROOOOOOOOOO!");
			return data[D*H*i + D*j + k];
		}
		void set(int i, int j, int k, double val) {
			if (i > W)
				printf("ERRRRRROOOOOOOOOO!");
			if (j > H)
				printf("ERRRRRROOOOOOOOOO!");
			if (k > D)
				printf("ERRRRRROOOOOOOOOO!");
			data[D*H*i + D*j + k] = val;
		}
		void add(int i, int j, int k, double val) {
			if (i > W)
				printf("ERRRRRROOOOOOOOOO!");
			if (j > H)
				printf("ERRRRRROOOOOOOOOO!");
			if (k > D)
				printf("ERRRRRROOOOOOOOOO!");
			data[D*H*i + D*j + k] += val;
		}
	};

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
				W_count = 1;
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
		void fit(SequenceVector &data_seq,  SequenceVector &gen_seq,int max_iter,bool useConjugate);
		void evalfP_gen(SequenceVector &gen_seq);
		double evalfMaxDelta(SequenceVector &gen_seq,int Lsize,int VJsize,int AAsize,vector<double>& data_Ldistibution, Martrix2d& data_VJdistibution, Martrix3d& data_AAdistibution);
		double evalfMeanDelta(SequenceVector &gen_seq, int Lsize, int VJsize, int AAsize, vector<double>& data_Ldistibution, Martrix2d& data_VJdistibution, Martrix3d& data_AAdistibution);
		double evalfMaxLikehood(SequenceVector &gen_seq,SequenceVector &data_seq,double alpha, vector<double>& F,int Lsize,int VJsize,int AAsize);
		double optimizeStep(SequenceVector &gen_seq,SequenceVector &data_seq, vector<double>& F, vector<double>& emptyArray,int Lsize,int VJsize,int AAsize,double* maxDelt);

        /**
         *
         */
		double predict(const Sequence &seq);
		double* predictMany(const SequenceVector &seq);

		double* Gw ;
		double* Dw ;

		vector<double>& get_q_L(){return q_L;};
		vector<double>& get_q_VJ(){return q_VJ.data;};
		vector<double>& get_ilA(){return q_ilA.data;};


		vector<double>& get_Pgen_L(){return gen_L;};
		vector<double>& get_Pgen_VJ(){	return gen_VJ.data;};
		vector<double>& get_Pgen_ilA(){return gen_ilA.data;};
		double getZ(){return Z;};
		int getMinL(){return minL;};
		int getMaxL(){return maxL;};

		map<string,int>* getVIndexes(){return V_indexes;};
		map<string,int>* getJIndexes(){return J_indexes;};
		map<char,int>* getAAindexes(){return &aminoAcidIndexes;};

        
	private:
		static const double EPS;
		static const int MAX_STEP;
		static const int MinCountForLength;
		static const int ScalarOptimizeN;

		int minL,maxL;
		map<string,int>* V_indexes, * J_indexes;
		map<char,int> aminoAcidIndexes;



		vector<double> q_L;
		Martrix2d q_VJ;
		Martrix3d q_ilA;

		
		

		vector<double> gen_L;
		Martrix2d gen_VJ;
		Martrix3d gen_ilA;



		double Z;

		double Q(const Sequence &seq);

	
		void findMinMaxLength(const SequenceVector &data_seq, const SequenceVector &gen_seq);
		void evalfDataLDistribution(const SequenceVector &data_seq, int minL, int maxL, int minFrequency, vector<double>& l_distribution);
		map<string,int>* extractVSet(const SequenceVector &data_seq, const SequenceVector &gen_seq);
		map<string,int>* extractJSet(const SequenceVector &data_seq, const SequenceVector &gen_seq);

		void evalf_gen_Ldistribution(const SequenceVector &gen_seq, vector<double>& l_distribution);
		double evalf_Z(const SequenceVector &gen_seq);

		void transformData(SequenceVector *seq);
	};
}

