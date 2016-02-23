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
		internal.push_back(tok);
	}
  
}

void parse_data_file(char* path,vector<Sequence>& out){
	out.reserve(200000);
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
	out.reserve(200000);
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

int _tmain(int argc, _TCHAR* argv[])
{
	vector<Sequence> gen;
	vector<Sequence> data;
	
	cout<<"start parsing data file\n";
	parse_data_file("C://immunology//github//mimir//build//Debug//data.txt",data);
	cout<<"start parsing gen file\n";
	parse_gen_file("C://immunology//github//mimir//build//Debug//gen.txt",gen);

	SelectionModel* S=new SelectionModel();
	cout<<"start fiting\n";
	S->fit(data,gen);
	cout<<"fit down";

	return 0;
}

