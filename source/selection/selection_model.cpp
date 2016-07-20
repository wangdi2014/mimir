#include "selection_model.h"
#include<time.h>


using namespace mimir;

const double SelectionModel::EPS=0.001;
const int SelectionModel::MAX_STEP=300;

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

void SelectionModel::evalfP_gen(SequenceVector &gen_seq){
	const int Lsize = maxL - minL + 1;
	const int VJsize = V_indexes->size()*J_indexes->size();
	const int AAsize = (maxL - minL + 1)*maxL*aminoAcidIndexes.size();
	const int fullSize=Lsize+VJsize+AAsize;

	gen_L=new double[maxL-minL+1];
	memset(gen_L,0,sizeof(double)*(maxL-minL+1));
	gen_VJ=new double[V_indexes->size()*J_indexes->size()];
	memset(gen_VJ,0,sizeof(double)*(V_indexes->size()*J_indexes->size()));
	gen_ilA=new double[(maxL-minL+1)*maxL*aminoAcidIndexes.size()];
	memset(gen_ilA,0,sizeof(double)*((maxL-minL+1)*maxL*aminoAcidIndexes.size()));

	float Qsum=0;
	for (int i = 0; i < gen_seq.size(); i++){

			int v_in = gen_seq[i].Vindex[0];
			int j_in = gen_seq[i].Jindex[0];
			if(gen_seq[i].aminoacide.length()<minL||gen_seq[i].aminoacide.length()>maxL)
				continue;
			int Lindex = gen_seq[i].aminoacide.length() - minL;
			int VJindex=v_in*J_indexes->size() + j_in;
			int L = gen_seq[i].aminoacide.length();
			
			Qsum+=1;
			gen_L[Lindex] += 1;
			gen_VJ[VJindex] += 1;
			for (int j = 0; j<L; j++){
				int aa_index = gen_seq[i].aminoacide_indexes[j] * (maxL - minL + 1)*maxL + (L - minL)*maxL + j;
				gen_ilA[aa_index] += 1;
			}
	}
	printf("Qsum::::%f",Qsum);
	for(int i=0;i<Lsize;i++){
		gen_L[i]=gen_L[i]/Qsum;
	}
	for(int i=0;i<VJsize;i++){
		gen_VJ[i]=gen_VJ[i]/Qsum;
	}
	for(int i=0;i<AAsize;i++){
		gen_ilA[i]=gen_ilA[i]/Qsum;
	}
}

void SelectionModel::fit(SequenceVector &data_seq, SequenceVector &gen_seq,int max_iter){
	srand(time(NULL));
	findMinMaxLength(data_seq, gen_seq);
	
	data_Ldistribution = evalfDataLDistribution(data_seq, minL, maxL, 1000);
	int newMinL,newMaxL;
	for(int l=minL;l<=maxL;l++){
		if(data_Ldistribution[l-minL]>1000){
			newMinL=l;
			break;
		}
	}
	for(int l=maxL;l>=minL;l--){
		if(data_Ldistribution[l-minL]>1000){
			newMaxL=l;
			break;
		}
	}
	printf("minL %d,maxL %d",newMinL,newMaxL);
	minL=newMinL;
	maxL=newMaxL;
	delete[] data_Ldistribution;
	data_Ldistribution=new double[maxL-minL+1];
	memset(data_Ldistribution,0,sizeof(double)*(maxL-minL+1));
	int N = 0;
	long total_sum = 0;
	for(int i=0;i<data_seq.size();i++){
		if(data_seq[i].aminoacide.size()>=minL&&data_seq[i].aminoacide.size()<=maxL){
			N+=1;
			total_sum += data_seq[i].W_count;
			data_Ldistribution[data_seq[i].aminoacide.size()-minL]+=data_seq[i].W_count;
		}
	}
	

	V_indexes = extractVSet(data_seq, gen_seq);
	J_indexes = extractJSet(data_seq, gen_seq);

	
	transformData(&gen_seq);
	transformData(&data_seq);



	data_VJpairDistribution = new double[V_indexes->size()*J_indexes->size()];
	memset(data_VJpairDistribution, 0, sizeof(double)*V_indexes->size()*J_indexes->size());


	for (int i = 0; i < data_seq.size(); i++){
		for(int v_counter=0;v_counter<data_seq[i].V.size();v_counter++){
			for(int j_counter=0;j_counter<data_seq[i].J.size();j_counter++){
				int v_in = data_seq[i].Vindex[v_counter];
				int j_in = data_seq[i].Jindex[j_counter];

				float w = data_seq[i].weights[v_counter*data_seq[i].J.size() + j_counter] * data_seq[i].W_count;
				
				if(data_seq[i].aminoacide.length()>=minL&&data_seq[i].aminoacide.length()<=maxL){
					data_VJpairDistribution[v_in*J_indexes->size() + j_in] += w;
				}
			}
		}
	}
	for (int i = 0; i < V_indexes->size()*J_indexes->size(); i++)
		data_VJpairDistribution[i] /= (double)total_sum;


	//initialize aaDistribution
	data_AAdistibution = new double[(maxL - minL + 1)*maxL*aminoAcidIndexes.size()];
	memset(data_AAdistibution, 0, sizeof(double)*(maxL - minL + 1)*maxL*aminoAcidIndexes.size());
	for (int i = 0; i < data_seq.size(); i++){
		int L = data_seq[i].aminoacide.length();
		if(L>=minL&&L<=maxL){
			for (int j = 0; j < L; j++){
				char sym = data_seq[i].aminoacide[j];
				if (sym != '*'&&sym != '~'){
					int aa_index = aminoAcidIndexes.at(sym);
					data_AAdistibution[aa_index*(maxL - minL + 1)*maxL + (L - minL)*maxL + j] += data_seq[i].W_count;
				}
			}
		}
	}
	for (int i = 0; i < maxL - minL + 1; i++) 
		for(int j=0;j<maxL;j++)
			for(int k=0;k<aminoAcidIndexes.size();k++)
				data_AAdistibution[k*(maxL - minL + 1)*maxL+i*maxL+j] /= data_Ldistribution[i];
	

	for (int i = 0;i<newMaxL - newMinL + 1;i++) {
		data_Ldistribution[i] /= (double)total_sum;
	}

	const int Lsize = maxL - minL + 1;
	const int VJsize = V_indexes->size()*J_indexes->size();
	const int AAsize = (maxL - minL + 1)*maxL*aminoAcidIndexes.size();
	const int fullSize=Lsize+VJsize+AAsize;

	//initialize start parametr
	q_L=new double[maxL-minL+1];
	for(int i=0;i<(maxL-minL+1);i++) q_L[i]=((double) rand() / (RAND_MAX))+0.5;
	//memset(q_L,0,sizeof(double)*(maxL-minL+1));
	q_VJ=new double[V_indexes->size()*J_indexes->size()];
	//memset(q_VJ,0,sizeof(double)*V_indexes->size()*J_indexes->size());
	for(int i=0;i<(V_indexes->size()*J_indexes->size());i++) q_VJ[i]=((double) rand() / (RAND_MAX))+0.5;
	q_ilA=new double[(maxL-minL+1)*maxL*aminoAcidIndexes.size()];
	//memset(q_ilA,0,sizeof(double)*(maxL-minL+1)*(maxL-minL+1)*aminoAcidIndexes.size());
	for(int i=0;i<((maxL-minL+1)*maxL*aminoAcidIndexes.size());i++) q_ilA[i]=((double) rand() / (RAND_MAX))+0.5;
	Z=1.0;


	//initialize gen distribution
	double* F=new double[fullSize];
	double* emptyArray=new double[fullSize];

	


	//start iteration
	double delt=1;
	int counter=0;
	
	//main loop
	double gradient_step=1;



	

	evalfP_gen(gen_seq);
	ofstream file;
	file.open("C://immunology//github//mimir//build//Debug//log_likehood_4");

	while(counter<max_iter){
		counter++;
		double newVal;
		delt=0;
		double Qsum=0;
		memset(F,0,fullSize*sizeof(double));
		for (int i = 0; i < gen_seq.size(); i++){

			int v_in = gen_seq[i].Vindex[0];
			int j_in = gen_seq[i].Jindex[0];
			if(gen_seq[i].aminoacide.length()<minL||gen_seq[i].aminoacide.length()>maxL)
				continue;
			int Lindex = gen_seq[i].aminoacide.length() - minL;
			int VJindex=Lsize+v_in*J_indexes->size() + j_in;
			int L = gen_seq[i].aminoacide.length();
			double q = 1;

			for (int k = 0; k < L; k++){
				int aa_index = gen_seq[i].aminoacide_indexes[k];
				q *= q_ilA[aa_index*(maxL - minL + 1)*maxL + Lindex*maxL + k];

			}
			q*=q_VJ[v_in*J_indexes->size() + j_in] * q_L[gen_seq[i].aminoacide.length() - minL];
			Qsum+=q;
			F[Lindex] += q;
			F[VJindex] += q;
			for (int j = 0; j<L; j++){
				int aa_index = Lsize+VJsize+gen_seq[i].aminoacide_indexes[j] * (maxL - minL + 1)*maxL + (L - minL)*maxL + j;
				F[aa_index] += q;
			}
		}
		
		for (int i = 0; i < maxL - minL + 1; i++)
			for (int j = 0;j<maxL;j++)
				for (int k = 0;k < aminoAcidIndexes.size();k++)
				{
					int index = k*(maxL - minL + 1)*maxL + i*maxL + j+Lsize+VJsize;
					F[index] /= F[i];
				}

		for (int i = 0;i<Lsize + VJsize;i++) {
			F[i] /= Qsum;
		}

		for(int i=0;i<fullSize;i++){
			
			if(F[i]!=0){
				if(i<Lsize){
					F[i]=(data_Ldistribution[i]-F[i])/q_L[i];
					continue;
				}
				if(i<Lsize+VJsize){
					F[i]=(data_VJpairDistribution[i-Lsize]-F[i])/q_VJ[i-Lsize];
					continue;
				}
				F[i]=(data_AAdistibution[i-Lsize-VJsize]-F[i])/q_ilA[i-Lsize-VJsize];
			}
			else{
			}
		}
 		gradient_step=optimizeStep(gen_seq,data_seq,F,emptyArray,Lsize,VJsize,AAsize,&delt);
		
		printf("step:::%f",gradient_step);
		for(int i=0;i<fullSize;i++){
			if(F[i]!=0){
				if(i<Lsize){
					q_L[i]+=gradient_step*F[i];
					if(q_L[i]<0) q_L[i]=0;
					continue;
				}
				if(i<Lsize+VJsize){
					q_VJ[i-Lsize]+=gradient_step*F[i];
					if(q_VJ[i-Lsize]<0) q_VJ[i-Lsize]=0;

					continue;
				}
				q_ilA[i-Lsize-VJsize]+=gradient_step*F[i];
				if(q_ilA[i-Lsize-VJsize]<0) q_ilA[i-Lsize-VJsize]=0;
			}
		}
		Z=evalf_Z(gen_seq);
		file<<to_string(delt)+"\n";
		printf("delt %f \n",delt);
	}

	int i;
	for(i=0;i<Lsize;i++){
		if(gen_L[i]==0)
			q_L[i]=0;
	}
	for(i=0;i<VJsize;i++){
		if(gen_VJ[i]==0)
			q_VJ[i]=0;
	}
	for(i=0;i<AAsize;i++){
		if(gen_ilA[i]==0)
			q_ilA[i]=0;
	}
	file.close();
	printf("delt %f \n",delt);
	printf("step %d ",counter);


}

double SelectionModel::evalfMaxDelta(SequenceVector &gen_seq,double alpha,double* F,double* emptyArray,int Lsize,int VJsize,int AAsize){
	int fullSize=Lsize+VJsize+AAsize;
	double Qsum=0;
	memset(emptyArray,0,sizeof(double)*fullSize);
	for(int i=0;i<gen_seq.size();i++){
		int v_in = gen_seq[i].Vindex[0];
			int j_in = gen_seq[i].Jindex[0];
			if(gen_seq[i].aminoacide.length()<minL||gen_seq[i].aminoacide.length()>maxL)
				continue;
			int Lindex = gen_seq[i].aminoacide.length() - minL;
			int VJindex=Lsize+v_in*J_indexes->size() + j_in;
			int L = gen_seq[i].aminoacide.length();
			double q = 1;

			for (int k = 0; k < L; k++){
				int aa_index = gen_seq[i].aminoacide_indexes[k];
				q *= q_ilA[aa_index*(maxL - minL + 1)*maxL + Lindex*maxL + k]+alpha*F[Lsize+VJsize+aa_index*(maxL - minL + 1)*maxL + Lindex*maxL + k];

			}
			q*=(q_VJ[v_in*J_indexes->size() + j_in]+alpha*F[Lsize+v_in*J_indexes->size() + j_in]) * (q_L[gen_seq[i].aminoacide.length() - minL]+alpha*F[gen_seq[i].aminoacide.length() - minL]);
			Qsum+=q;
			emptyArray[Lindex] += q;
			emptyArray[VJindex] += q;
			for (int j = 0; j<L; j++){
				int aa_index = Lsize+VJsize+gen_seq[i].aminoacide_indexes[j] * (maxL - minL + 1)*maxL + (L - minL)*maxL + j;
				emptyArray[aa_index] += q;
			}
	}
	double maxDelt=0;
	int i=0;
	for(i=0;i<Lsize;i++){
		if(emptyArray[i]!=0){
			double d=abs(emptyArray[i]/Qsum-data_Ldistribution[i])/q_L[i];
			if(d>maxDelt) maxDelt=d;
		}
	}
	for(;i<VJsize+Lsize;i++){
		if(emptyArray[i]!=0){
			double d=abs(emptyArray[i]/Qsum-data_VJpairDistribution[i-Lsize])/q_VJ[i-Lsize];
		if(d>maxDelt) maxDelt=d;
		}
	}
	for(;i<fullSize;i++){
		if(emptyArray[i]!=0){
			double d=abs(emptyArray[i]/Qsum-data_AAdistibution[i-Lsize-VJsize])/q_ilA[i-Lsize-VJsize];
		if(d>maxDelt) maxDelt=d;
		}
	}

	return maxDelt;
}

double SelectionModel::evalfMaxLikehood(SequenceVector &gen_seq,SequenceVector &data_seq,double alpha,double* F,int Lsize,int VJsize,int AAsize){
	int fullSize=Lsize+VJsize+AAsize;
	double Qsum=0;
	for(int i=0;i<gen_seq.size();i++){
				int v_in = gen_seq[i].Vindex[0];
				int j_in = gen_seq[i].Jindex[0];
				if(gen_seq[i].aminoacide.length()<minL||gen_seq[i].aminoacide.length()>maxL)
					continue;
				int Lindex = gen_seq[i].aminoacide.length() - minL;
				int VJindex=Lsize+v_in*J_indexes->size() + j_in;
				int L = gen_seq[i].aminoacide.length();
				double q = 1;

				for (int k = 0; k < L; k++){
					int aa_index = gen_seq[i].aminoacide_indexes[k];
					q *= q_ilA[aa_index*(maxL - minL + 1)*maxL + Lindex*maxL + k]+alpha*F[Lsize+VJsize+aa_index*(maxL - minL + 1)*maxL + Lindex*maxL + k];

				}
				q*=(q_VJ[v_in*J_indexes->size() + j_in]+alpha*F[Lsize+v_in*J_indexes->size() + j_in]) * (q_L[gen_seq[i].aminoacide.length() - minL]+alpha*F[gen_seq[i].aminoacide.length() - minL]);
				Qsum+=q;
	}
	double Z=Qsum/gen_seq.size();
	double likehood=0;
	for(int counter=0;counter<data_seq.size();counter++){
		int i=counter;
		if(data_seq[i].valid){
			if(data_seq[i].aminoacide.length()<minL||data_seq[i].aminoacide.length()>maxL)
					continue;
			for(int v_counter=0;v_counter<data_seq[i].V.size();v_counter++){
				for(int j_counter=0;j_counter<data_seq[i].J.size();j_counter++){

						float w=data_seq[i].weights[v_counter*data_seq[i].J.size()+j_counter]*data_seq[i].W_count;
						int v_in = data_seq[i].Vindex[v_counter];
						int j_in = data_seq[i].Jindex[j_counter];
				
						int Lindex = data_seq[i].aminoacide.length() - minL;
						int VJindex=Lsize+v_in*J_indexes->size() + j_in;
						int L = data_seq[i].aminoacide.length();
						double q = 1;

						for (int k = 0; k < L; k++){
							int aa_index = data_seq[i].aminoacide_indexes[k];
							q *= q_ilA[aa_index*(maxL - minL + 1)*maxL + Lindex*maxL + k]+alpha*F[Lsize+VJsize+aa_index*(maxL - minL + 1)*maxL + Lindex*maxL + k];
						}
						q*=(q_VJ[v_in*J_indexes->size() + j_in]+alpha*F[Lsize+v_in*J_indexes->size() + j_in]) * (q_L[data_seq[i].aminoacide.length() - minL]+alpha*F[data_seq[i].aminoacide.length() - minL]);
						likehood+=w*(log(q)-log(Z));
				}
			}
		}
	}
	return likehood;

}

double SelectionModel::optimizeStep(SequenceVector &gen_seq,SequenceVector &data_seq,double* F,double* emptyArray,int Lsize,int VJsize,int AAsize,double* delt){
	double phi=0.5*(1+sqrt(5));
	double a=0;
	double b=10;
	while(abs(b-a)>0.01){
		double x1=b-(b-a)/phi;
		double x2=a+(b-a)/phi;
		double y1=evalfMaxLikehood(gen_seq,data_seq,x1,F,Lsize,VJsize,AAsize);
		double y2=evalfMaxLikehood(gen_seq,data_seq,x2,F,Lsize,VJsize,AAsize);
		if(y1<=y2){
			a=x1;
			x1=x2;
			x2=a+(b-a)/phi;
		}
		else
		{
			b=x2;
			x2=x1;
			x1=b-(b-a)/phi;
		}
	}
	double out=0.5*(a+b);
	//double out=1;
	*delt=evalfMaxLikehood(gen_seq,data_seq,out,F,Lsize,VJsize,AAsize);
	return out;
}

double SelectionModel::predict(const Sequence &seq){
	return Q(seq);
}

double* SelectionModel::predictMany(const SequenceVector &seq){
	double* prediction=new double(seq.size());
	for(int i=0;i<seq.size();i++){
		prediction[i]=Q(seq[i]);
	}
	return prediction;
}

double SelectionModel::Q(const Sequence &seq){
	int v=seq.Vindex[0];
	int j=seq.Jindex[0];
	int Lindex=seq.aminoacide.length()-minL;
	double q=q_VJ[v*J_indexes->size()+j]*q_L[seq.aminoacide.length()-minL];
	
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

double* SelectionModel::evalfDataLDistribution(const SequenceVector &data_seq,int minL,int maxL,int minFrequency){
	double* LDistribution=new double[maxL-minL+1];
	memset(LDistribution,0,sizeof(double)*(maxL-minL+1));
	int N=data_seq.size();
	for(int i=0;i<N;i++){
		int L=data_seq[i].aminoacide.length();
		LDistribution[L-minL]+=1;
	}
	return LDistribution;
}

map<string,int>* SelectionModel::extractVSet(const SequenceVector &data_seq,const SequenceVector &gen_seq){
	map<string,int>* out=new map<string,int>();
	int N=data_seq.size();
	int k=0;
	for(int i=0;i<N;i++){
		for(int v_counter=0;v_counter<data_seq[i].V.size();v_counter++){
			if(out->find(data_seq[i].V[v_counter]) == out->end() ){
				(*out)[data_seq[i].V[v_counter]]=k;
				k++;
			}
		}
	}
	N=gen_seq.size();
	for(int i=0;i<N;i++){
		if(out->find(gen_seq[i].V[0]) == out->end() ){
			(*out)[gen_seq[i].V[0]]=k;
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
		for(int j_counter=0;j_counter<data_seq[i].J.size();j_counter++){
			if(out->find(data_seq[i].J[j_counter]) == out->end() ){
				(*out)[data_seq[i].J[j_counter]]=k;
				k++;
			}
		}
	}
	N=gen_seq.size();
	for(int i=0;i<N;i++){
		if(out->find(gen_seq[i].J[0]) == out->end() ){
			(*out)[gen_seq[i].J[0]]=k;
			k++;
		}
	}
	return out;
}

void SelectionModel::evalf_gen_Ldistribution(const SequenceVector &gen_seq,double* l_distribution){
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

void SelectionModel::evalf_gen_VJdistribution(const SequenceVector &gen_seq,double* VJ_distribution){
	memset(VJ_distribution,0,sizeof(double)*V_indexes->size()*J_indexes->size());

	double sum=0;
	for(int i=0;i<gen_seq.size();i++){
		if(gen_seq[i].valid){
			int vindex=gen_seq[i].Vindex[0];
			int jindex=gen_seq[i].Jindex[0];
			double q=Q(gen_seq[i]);
			VJ_distribution[vindex*J_indexes->size()+jindex]+=q;
			sum+=q;
		}
	}
	for(int i=0;i<V_indexes->size()*J_indexes->size();i++){
		VJ_distribution[i]/=sum;
	}
}

void SelectionModel::evalf_gen_AAdistribution(const SequenceVector &gen_seq,double* AA_distribution){
	memset(AA_distribution,0,sizeof(double)*aminoAcidIndexes.size()*maxL*(maxL-minL+1));

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

double SelectionModel::evalf_Z(const SequenceVector &gen_seq){
	double sum=0;
	for(int i=0;i<gen_seq.size();i++){
		if(gen_seq[i].valid){
			if(gen_seq[i].aminoacide.length()<minL||gen_seq[i].aminoacide.length()>maxL)
				continue;
			int v=gen_seq[i].Vindex[0];
			int j=gen_seq[i].Jindex[0];
			int Lindex=gen_seq[i].aminoacide.length()-minL;
			double q=q_VJ[v*J_indexes->size()+j]*q_L[gen_seq[i].aminoacide.length()-minL];
	
			for(int k=0;k<gen_seq[i].aminoacide.length();k++){
				int aa_index=gen_seq[i].aminoacide_indexes[k];
				q*=q_ilA[aa_index*(maxL-minL+1)*maxL+Lindex*maxL+k];
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