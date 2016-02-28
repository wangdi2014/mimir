#include "selection_model.h"


using namespace mimir;

const float SelectionModel::EPS=0.0000001;
const int SelectionModel::MAX_STEP=10;

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
	//aa_map['*']=20;
	return aa_map;
	
}




SelectionModel::SelectionModel(void)
{
	V_indexes=NULL;
	J_indexes=NULL;

	data_Ldistribution=NULL;
	q_L=NULL;

	data_VJpairDistribution=NULL;
	q_VJ=NULL;

	data_AAdistibution=NULL;
	q_ilA=NULL;

	minL=99999;
	maxL=0;

	SelectionModel::aminoAcidIndexes=init_aa_map();
}

SelectionModel::~SelectionModel(void)
{
	delete V_indexes;
	delete J_indexes;
	delete[] data_Ldistribution;
	delete[] q_L;
	
	delete[] data_VJpairDistribution;
	delete[] q_VJ;

	delete[] data_AAdistibution;
	delete[] q_ilA;
}

void SelectionModel::fit(const SequenceVector &data_seq, SequenceVector &gen_seq){
	findMinMaxLength(data_seq,gen_seq);
	data_Ldistribution=evalfDataLDistribution(data_seq,minL,maxL,1000);
	int N=data_seq.size();
	V_indexes=extractVSet(data_seq,gen_seq);
	J_indexes=extractJSet(data_seq,gen_seq);

	data_VJpairDistribution=new float[V_indexes->size()*J_indexes->size()];
	memset(data_VJpairDistribution,0,sizeof(float)*V_indexes->size()*J_indexes->size());


	for(int i=0;i<N;i++){
		int v_index=(*V_indexes)[data_seq[i].V];
		int j_index=(*J_indexes)[data_seq[i].J];

		data_VJpairDistribution[v_index*J_indexes->size()+j_index]+=1;
	}
	for(int i=0;i<V_indexes->size()*J_indexes->size();i++)
		data_VJpairDistribution[i]/=(float)N;
		
	//initialize aaDistribution
	data_AAdistibution=new float[(maxL-minL+1)*maxL*aminoAcidIndexes.size()];
	memset(data_AAdistibution,0,sizeof(float)*(maxL-minL+1)*maxL*aminoAcidIndexes.size());
	for(int i=0;i<N;i++){
		int L=data_seq[i].aminoacide.length();
		for(int j=0;j<L;j++){
			char sym=data_seq[i].aminoacide[j];
			if(sym!='*'&&sym!='~'){
				int aa_index=aminoAcidIndexes.at(sym);
				data_AAdistibution[aa_index*(maxL-minL+1)*maxL+(L-minL)*maxL+j]+=1;
			}
		}
	}
	for(int i=0;i<(maxL-minL+1)*maxL*aminoAcidIndexes.size();i++)
		data_AAdistibution[i]/=(float)N;

	//initialize start parametr
	q_L=new float[maxL-minL+1];
	for(int i=0;i<(maxL-minL+1);i++) q_L[i]=1;
	//memset(q_L,0,sizeof(float)*(maxL-minL+1));
	q_VJ=new float[V_indexes->size()*J_indexes->size()];
	//memset(q_VJ,0,sizeof(float)*V_indexes->size()*J_indexes->size());
	for(int i=0;i<(V_indexes->size()*J_indexes->size());i++) q_VJ[i]=1;
	q_ilA=new float[(maxL-minL+1)*maxL*aminoAcidIndexes.size()];
	//memset(q_ilA,0,sizeof(float)*(maxL-minL+1)*(maxL-minL+1)*aminoAcidIndexes.size());
	for(int i=0;i<((maxL-minL+1)*maxL*aminoAcidIndexes.size());i++) q_ilA[i]=1;
	Z=1.0;


	//initialize gen distribution
	float* gen_Ldistribution=new float[maxL-minL+1];
	float* gen_VJdistribution=new float[V_indexes->size()*J_indexes->size()];
	float* gen_AAdistribution=new float[(maxL-minL+1)*maxL*aminoAcidIndexes.size()];

	//start iteration
	float delt=1;
	int counter=0;
	
	//main loop
	float gradient_step=1;

	transformData(&gen_seq);

	while(delt>EPS&&counter<MAX_STEP){
		counter++;
		float newVal;
		delt=0;

		evalf_gen_Ldistribution(gen_seq,gen_Ldistribution);
		for(int i=0;i<(maxL-minL+1);i++){
			newVal=q_L[i]+gradient_step*(data_Ldistribution[i]-gen_Ldistribution[i]);
			if(abs(newVal-q_L[i])>delt)
				delt=abs(newVal-q_L[i]);
			q_L[i]=newVal;
		}
		
		evalf_gen_VJdistribution(gen_seq,gen_VJdistribution);
		for(int i=0;i<V_indexes->size()*J_indexes->size();i++){
			if(gen_VJdistribution[i]<0){
				printf("1");
			}
			newVal=q_VJ[i]+gradient_step*(data_VJpairDistribution[i]-gen_VJdistribution[i]);
			if(abs(newVal-q_VJ[i])>delt)
				delt=abs(newVal-q_VJ[i]);
			q_VJ[i]=newVal;
		}
		evalf_gen_AAdistribution(gen_seq,gen_AAdistribution);
		for(int i=0;i<(maxL-minL+1)*maxL*aminoAcidIndexes.size();i++){
			
			newVal=q_ilA[i]+gradient_step*(data_AAdistibution[i]-gen_AAdistribution[i]);
			if(abs(newVal-q_ilA[i])>delt)
				delt=abs(newVal-q_ilA[i]);
			q_ilA[i]=newVal;
		}
		Z=evalf_Z(gen_seq);
		printf("delt %f \n",delt);
	}
	printf("delt %f \n",delt);
	printf("step %d ",counter);
	delete[] gen_Ldistribution;
	delete[] gen_VJdistribution;
	delete[] gen_AAdistribution;

}

float SelectionModel::predict(const Sequence &seq){
	return Q(seq);
}

float* SelectionModel::predictMany(const SequenceVector &seq){
	float* prediction=new float(seq.size());
	for(int i=0;i<seq.size();i++){
		prediction[i]=Q(seq[i]);
	}
	return prediction;
}








float SelectionModel::Q(const Sequence &seq){
	int v=seq.Vindex;
	int j=seq.Jindex;
	int Lindex=seq.aminoacide.length()-minL;
	float q=q_VJ[v*J_indexes->size()+j]*q_L[seq.aminoacide.length()-minL];
	
	for(int i=0;i<seq.aminoacide.length();i++){
		int aa_index=seq.aminoacide_indexes[i];
		q*=q_ilA[aa_index*(maxL-minL+1)*(maxL-minL+1)+Lindex*maxL+i];
	}
	return q/Z;
}

void SelectionModel::findMinMaxLength(const SequenceVector &data_seq,const SequenceVector &gen_seq){
	int N=data_seq.size();
	for(int i=0;i<N;i++){
		int L=data_seq[i].aminoacide.length();
		if(L>maxL) maxL=L;
		if(L<minL) minL=L;
	}
	N=gen_seq.size();
	for(int i=0;i<N;i++){
		int L=gen_seq[i].aminoacide.length();
		if(L>maxL) maxL=L;
		if(L<minL) minL=L;
	}
}

float* SelectionModel::evalfDataLDistribution(const SequenceVector &data_seq,int minL,int maxL,int minFrequency){
	float* LDistribution=new float[maxL-minL+1];
	memset(LDistribution,0,sizeof(float)*(maxL-minL+1));
	int N=data_seq.size();
	for(int i=0;i<N;i++){
		int L=data_seq[i].aminoacide.length();
		LDistribution[L-minL]+=1;
	}
	for(int i=0;i<maxL-minL+1;i++)
		if(LDistribution[i]>=minFrequency)
			LDistribution[i]/=(float)N;
		else
			LDistribution[i]=0;
	return LDistribution;
}

map<string,int>* SelectionModel::extractVSet(const SequenceVector &data_seq,const SequenceVector &gen_seq){
	map<string,int>* out=new map<string,int>();
	int N=data_seq.size();
	int k=0;
	for(int i=0;i<N;i++){
		if(out->find(data_seq[i].V) == out->end() ){
			(*out)[data_seq[i].V]=k;
			k++;
		}
	}
	N=gen_seq.size();
	for(int i=0;i<N;i++){
		if(out->find(gen_seq[i].V) == out->end() ){
			(*out)[gen_seq[i].V]=k;
			k++;
		}
	}
	return out;
}

map<string,int>* SelectionModel::extractJSet(const SequenceVector &data_seq,const SequenceVector &gen_seq){
	map<string,int>* out=new map<string,int>();
	int N=data_seq.size();
	int k=0;
	for(int i=0;i<N;i++){
		if(out->find(data_seq[i].J) == out->end() ){
			(*out)[data_seq[i].J]=k;
			k++;
		}
	}
	N=gen_seq.size();
	for(int i=0;i<N;i++){
		if(out->find(gen_seq[i].J) == out->end() ){
			(*out)[gen_seq[i].J]=k;
			k++;
		}
	}
	return out;
}

void SelectionModel::evalf_gen_Ldistribution(const SequenceVector &gen_seq,float* l_distribution){
	double* out=new double[maxL-minL+1];
	memset(out,0,sizeof(double)*(maxL-minL+1));
	double sum=0;
	for(int i=0;i<gen_seq.size();i++){
		if(gen_seq[i].valid){
			int l=gen_seq[i].aminoacide.size()-minL;
			double q=Q(gen_seq[i]);
			out[l]+=q;
			sum+=q;
		}
	}
	for(int i=0;i<maxL-minL+1;i++){
		l_distribution[i]=out[i]/sum;
	}
	delete[] out;
}

void SelectionModel::evalf_gen_VJdistribution(const SequenceVector &gen_seq,float* VJ_distribution){
	memset(VJ_distribution,0,sizeof(float)*V_indexes->size()*J_indexes->size());

	double sum=0;
	for(int i=0;i<gen_seq.size();i++){
		if(gen_seq[i].valid){
			int vindex=gen_seq[i].Vindex;
			int jindex=gen_seq[i].Jindex;
			double q=Q(gen_seq[i]);
			VJ_distribution[vindex*J_indexes->size()+jindex]+=q;
			sum+=q;
		}
	}
	for(int i=0;i<V_indexes->size()*J_indexes->size();i++){
		VJ_distribution[i]/=sum;
	}
}

void SelectionModel::evalf_gen_AAdistribution(const SequenceVector &gen_seq,float* AA_distribution){
	memset(AA_distribution,0,sizeof(float)*aminoAcidIndexes.size()*maxL*(maxL-minL+1));

	double sum=0;
	for(int i=0;i<gen_seq.size();i++){
		if(gen_seq[i].valid){
			double q=Q(gen_seq[i]);
			int L=gen_seq[i].aminoacide.size();
			for(int j=0;j<L;j++){
				AA_distribution[gen_seq[i].aminoacide_indexes[j]*(maxL-minL+1)*maxL+(L-minL)*maxL+j]+=q;
			}
			sum+=q;
		}
	}
	for(int i=0;i<aminoAcidIndexes.size()*maxL*(maxL-minL+1);i++){
		AA_distribution[i]/=sum;
	}
}

float SelectionModel::evalf_Z(const SequenceVector &gen_seq){
	double sum=0;
	for(int i=0;i<gen_seq.size();i++){
		if(gen_seq[i].valid){
			int v=gen_seq[i].Vindex;
			int j=gen_seq[i].Jindex;
			int Lindex=gen_seq[i].aminoacide.length()-minL;
			float q=q_VJ[v*J_indexes->size()+j]*q_L[gen_seq[i].aminoacide.length()-minL];
	
			for(int i=0;i<gen_seq[i].aminoacide.length();i++){
				int aa_index=gen_seq[i].aminoacide_indexes[i];
				q*=q_ilA[aa_index*(maxL-minL+1)*maxL+Lindex*maxL+i];
			}
			sum+=q;
		}
	}
	return sum/gen_seq.size();
}

void SelectionModel::transformData(SequenceVector *gen){
	for(int i=0;i<gen->size();i++){
		(*gen)[i].convertToIndexes(*V_indexes,*J_indexes,aminoAcidIndexes);
	}
}