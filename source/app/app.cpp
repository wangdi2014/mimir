// test.cpp: определяет точку входа для консольного приложения.
//
#include <fstream>
#include <string>
#include <sstream>


#include <stdio.h>
#include <tchar.h>
#include <iostream>
#include "selection_model.h"


using namespace std;
using namespace mimir;



void split(string &str, char delimiter,vector<string> &internal) {
	

	stringstream ss(str); // Turn the string into a stream.
	string tok;
  
	while(std::getline(ss, tok, delimiter)) {
		internal.push_back(tok.substr(0, tok.find(",", 0)));
	}
  
}

void parse_data_file(char* path,vector<Sequence>& out){
	out.reserve(350000);
	ifstream file(path);

    string line; 
	vector<string> words;
	
	if (file.is_open())
	{
		getline (file, line);
		while (!file.eof())
		{
			getline (file, line);
			split(line,'\t',words);
			if(words.size()>9)
				out.emplace_back(words[5],words[7],words[9]);
			words.clear();
		}
	}
}

void parse_gen_file(char* path,vector<Sequence>& out){
	out.reserve(350000);
	ifstream file(path);
    string line; 
	vector<string> words;
	
	if (file.is_open())
	{
		getline (file, line);
		while (!file.eof())
		{
			getline (file, line);
			split(line,'\t',words);
			if(words.size()>2)
				out.emplace_back(words[0],words[1],words[2]);
			words.clear();
			
		}
	}
}

void writeToSingleJson(SelectionModel &model,char* path){
	ofstream file;
	double* L=model.get_q_L();
	file.open(path);
	file<<"{\n\t\"L\":{\n\t\t";
	int minL=model.getMinL();
	int maxL=model.getMaxL();
	file<<"\"minL\":"+to_string(minL)+",\n\t\t\"maxL\":"+to_string(maxL)+",\n\t\t\"qL\":["+to_string(L[0]);
	for(int i=1;i<maxL-minL+1;i++){
		file<<','+to_string(L[i]);
	}
	file<<"]";
	file<<"\n\t},";
	//write VJ
	map<string,int>* V_indexes=model.getVIndexes();

	file<<"\n\t\"VJ\":{\n";
	file<<"\t\t\"V_indexes\":{\""+V_indexes->begin()->first+"\":"+to_string(V_indexes->begin()->second);
	bool first=true;
	for(map<string,int>::iterator it=V_indexes->begin();it!=V_indexes->end();it++){
		if(!first){
			file<<",\""+it->first+"\":"+to_string(it->second);
		}
		else
			first=false;
	}
	file<<"},\n";

	map<string,int>* J_indexes=model.getJIndexes();
	file<<"\t\t\"J_indexes\":{\""+J_indexes->begin()->first+"\":"+to_string(J_indexes->begin()->second);
	first=true;
	for(map<string,int>::iterator it=J_indexes->begin();it!=J_indexes->end();it++){
		if(!first){
			file<<",\""+it->first+"\":"+to_string(it->second);
		}
		else
			first=false;
	}
	file<<"},\n";

	double* q_VJ=model.get_q_VJ();
	file<<"\t\t\"q_VJ\":[";

	first=true;
	for(int v_index=0;v_index<V_indexes->size();v_index++){
		if(!first)
			file<<",";
		else
			first=false;
		file<<"["+to_string(q_VJ[v_index*J_indexes->size()]);
		for(int j_index=1;j_index<J_indexes->size();j_index++){
			file<<","+to_string(q_VJ[v_index*J_indexes->size()+j_index]);
		}
		file<<"]";
	}
	file<<"]\n\t},\n";

	//write aminoacide
	file<<"\t\"q_Li\":{\n";
	map<char,int>* aminoAcidIndexes=model.getAAindexes();
	double* q_ilA=model.get_ilA();
	first=true;

	for(map<char,int>::iterator it=aminoAcidIndexes->begin();it!=aminoAcidIndexes->end();it++){
		if(!first)
			file<<",\n";
		else
			first=false;
		string aminoacide("a");
		aminoacide[0]=it->first;
		file<<"\t\t\""+aminoacide+"\":[";
		int a_index=it->second;
		bool firstline=true;
		for(int i=0;i<maxL-minL+1;i++){
			if(!firstline)
				file<<",";
			else
				firstline=false;
			file<<"["+to_string(q_ilA[a_index*(maxL-minL+1)*maxL+i*maxL]);
			for(int j=1;j<maxL;j++){
				file<<","+to_string(q_ilA[a_index*(maxL-minL+1)*maxL+i*maxL+j]);
			}
			file<<"]";
		}
		file<<"]";
	}
	file<<"\n\t},\n";
	//write Z
	file<<"\t\"Z\":"+to_string(model.getZ());

	file<<"\n}";


	file.close();
}

int _tmain(int argc, char* argv[])
{
	vector<Sequence> gen;
	vector<Sequence> data;
	
	cout<<"start parsing data file\n";
	parse_data_file(argv[1],data);
	cout<<"start parsing gen file\n";
	parse_gen_file(argv[2],gen);

	SelectionModel* S=new SelectionModel();
	cout<<"start fiting\n";
	S->fit(data,gen);
	cout<<"fit done";
	writeToSingleJson(*S,"coef.json");


	return 0;
}

