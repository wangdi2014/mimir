#include "selection_model.h"


using namespace mimir;

const float SelectionModel::EPS=0.1;
const int SelectionModel::MAX_STEP=10000;

map<char,int> init_aa_map(){
	map<char,int> aa_map;
	aa_map['G']=0;
	aa_map['A']=1;
	aa_map['V']=2;
	aa_map['I']=3;
	aa_map['L']=4;
	aa_map['P']=5;
	aa_map['S']=6;
	aa_map['T']=7;
	aa_map['C']=8;
	aa_map['M']=9;
	aa_map['D']=10;
	aa_map['N']=11;
	aa_map['E']=12;
	aa_map['Q']=13;
	aa_map['K']=14;
	aa_map['R']=15;
	aa_map['H']=16;
	aa_map['F']=17;
	aa_map['Y']=18;
	aa_map['W']=19;
	return aa_map;
	
}
map<char,int> aminoAcidIndexes=init_aa_map();



SelectionModel::SelectionModel(void)
{
	

	data_Ldistribution=NULL;
	q_L=NULL;

	data_VJpairDistribution=NULL;
	q_VJ=NULL;

	data_AAdistibution=NULL;
	q_ilA=NULL;

	minL=99999;
	maxL=0;

}

SelectionModel::~SelectionModel(void)
{
	delete[] data_Ldistribution;
	delete[] q_L;
	
	delete[] data_VJpairDistribution;
	delete[] q_VJ;

	delete[] data_AAdistibution;
	delete[] q_ilA;
}

void SelectionModel::fit(const SequenceVector &data_seq,const SequenceVector &gen_seq){
	findMinMaxLength(data_seq,gen_seq);
	data_Ldistribution=evalfDataLDistribution(data_seq,minL,maxL,1000);
	q_L=new float[maxL-minL+1];
	memset(q_L,0,sizeof(float)*(maxL-minL+1));

	V_indexes=extractVSet(data_seq,gen_seq);
	J_indexes=extractJSet(data_seq,gen_seq);

	data_VJpairDistribution=new float[V_indexes->size()*J_indexes->size()];
	memset(data_VJpairDistribution,0,sizeof(float)*V_indexes->size()*J_indexes->size());

	int N=data_seq.size();
	for(int i=0;i<N;i++){
		int v_index=(*V_indexes)[*data_seq[i].V];
		int j_index=(*J_indexes)[*data_seq[i].J];
		data_VJpairDistribution[v_index*J_indexes->size()+j_index]+=1;
	}
	for(int i=0;i<V_indexes->size()*J_indexes->size();i++)
		data_VJpairDistribution[i]/=(float)N;

	//initialize aaDistribution
	data_AAdistibution=new float[(maxL-minL+1)*(maxL-minL+1)*aminoAcidIndexes.size()];
	memset(data_AAdistibution,0,sizeof(float)*(maxL-minL+1)*(maxL-minL+1)*aminoAcidIndexes.size());
	for(int i=0;i<N;i++){
		string* V=data_seq[i].V;
		int L=V->length();
		for(int j=0;j<L;j++){
			int aa_index=aminoAcidIndexes[(*V)[j]];
			data_AAdistibution[aa_index*(maxL-minL+1)*(maxL-minL+1)+L*(maxL-minL+1)+j]+=1;
		}
	}
	for(int i=0;i<(maxL-minL+1)*(maxL-minL+1)*aminoAcidIndexes.size();i++)
		data_AAdistibution[i]/=(float)N;

	//initialize start parametr
	q_L=new float[(maxL-minL+1)];
	memset(q_L,1,sizeof(float)*(maxL-minL+1));
	q_VJ=new float[V_indexes->size()*J_indexes->size()];
	memset(q_VJ,1,sizeof(float)*V_indexes->size()*J_indexes->size());
	q_ilA=new float[(maxL-minL+1)*(maxL-minL+1)*aminoAcidIndexes.size()];
	memset(q_ilA,1,sizeof(float)*(maxL-minL+1)*(maxL-minL+1)*aminoAcidIndexes.size());


	//initialize gen distribution
	float* gen_Ldistribution=new float[maxL-minL+1];

	//start iteration
	float delt=1;
	int counter=0;
	
	//main loop
	while(delt>EPS&&counter<MAX_STEP){
		evalf_gen_Ldistribution(gen_seq,1000,gen_Ldistribution);
		for(int i=0;i<(maxL-minL+1);i++){
			
		}
	
	}


}

float SelectionModel::Q(const Sequence &seq){
	int v=(*V_indexes)[*(seq.V)];
	int j=(*J_indexes)[*(seq.J)];
	int Lindex=seq.aminoacide->length()-minL;
	float q=q_VJ[v*J_indexes->size()+j]*q_L[seq.aminoacide->length()-minL];
	
	for(int i=0;i<seq.aminoacide->length();i++){
		int aa_index=aminoAcidIndexes[(*seq.aminoacide)[i]];
		q*=q_ilA[aa_index*(maxL-minL+1)*(maxL-minL+1)+Lindex*(maxL-minL+1)+i];
	}
	return q/Z;
}

void SelectionModel::findMinMaxLength(const SequenceVector &data_seq,const SequenceVector &gen_seq){
	int N=data_seq.size();
	for(int i=0;i<N;i++){
		int L=data_seq[i].aminoacide->length();
		if(L>maxL) maxL=L;
		if(L<minL) minL=L;
	}
	N=gen_seq.size();
	for(int i=0;i<N;i++){
		int L=gen_seq[i].aminoacide->length();
		if(L>maxL) maxL=L;
		if(L<minL) minL=L;
	}
}

float* SelectionModel::evalfDataLDistribution(const SequenceVector &data_seq,int minL,int maxL,int minFrequency){
	float* LDistribution=new float[maxL-minL+1];
	memset(LDistribution,0,sizeof(float)*(maxL-minL+1));
	int N=data_seq.size();
	for(int i=0;i<N;i++){
		int L=data_seq[i].aminoacide->length();
		LDistribution[L-minL]+=1;
	}
	for(int i=0;i<N;i++)
		if(LDistribution[i]>=minFrequency)
			LDistribution[i]/=(float)N;
		else
			LDistribution[i]=0;
	return LDistribution;
}

map<string,int>* SelectionModel::extractVSet(const SequenceVector &data_seq,const SequenceVector &gen_seq){
	map<string,int>* out=new map<string,int>();
	size_t N;
	int k=0;
	for(int i=0;i<N;i++){
		if(out->find(*data_seq[i].V) == out->end() ){
			(*out)[*data_seq[i].V]=k;
			k++;
		}
	}
	N=data_seq.size();
	for(int i=0;i<N;i++){
		if(out->find(*gen_seq[i].V) == out->end() ){
			(*out)[*gen_seq[i].V]=k;
			k++;
		}
	}
	return out;
}

map<string,int>* extractJSet(const SequenceVector &data_seq,const SequenceVector &gen_seq){
	map<string,int>* out=new map<string,int>();
	int N=data_seq.size();
	int k=0;
	for(int i=0;i<N;i++){
		if(out->find(*data_seq[i].J) == out->end() ){
			(*out)[*data_seq[i].J]=k;
			k++;
		}
	}
	N=data_seq.size();
	for(int i=0;i<N;i++){
		if(out->find(*gen_seq[i].J) == out->end() ){
			(*out)[*gen_seq[i].J]=k;
			k++;
		}
	}
	return out;
}

void SelectionModel::evalf_gen_Ldistribution(const SequenceVector &gen_seq,int minFrequency,float* l_distribution){
	double* out=new double[maxL-minL+1];
	double sum=0;
	for(int i=0;i<gen_seq.size();i++){
		int l=gen_seq[i].aminoacide->size()-minL;
		double q=Q(gen_seq[i]);
		out[l]+=q;
		sum+=q;
	}
	for(int i=0;i<maxL-minL+1;i++){
		if(out[i]>=minFrequency)
			l_distribution[i]=out[i]/sum;
		else
			l_distribution[i]=0;
	}
	delete[] out;

}