#include "count_min_sketch.h"
#include <vector>
#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include "strict_fstream.hpp"
#include "zstr.hpp"




//TODO multithread



using namespace std;



struct unitig{
    string sequence;
    vector<string>colors;
};



struct list_unitig{
    vector<unitig> list;
    uint64_t unitigs_length;
};



bool is_fasta(const string& str){
    if (str.find(".fa") != std::string::npos) {
        return true;
    }
    if (str.find(".fasta") != std::string::npos) {
        return true;
    }
    return false;
}



uint64_t insert_fasta_in_filter(const string& filename,CountMinSketch& cms){
    zstr::ifstream f(filename);
    string ref, useless;
    uint64_t seq_number(0);
    while (not f.eof()) {
        getline(f, useless);   // read a comment, useless
        getline(f, ref);		// read the ACGT sequence
        if(not ref.empty()){
            cms_add(&cms,ref.c_str());
            // cout<<"i add:   "<<ref<<endl;
        }
        ref=useless="";
        seq_number++;
    }
    return seq_number;
}



void add_index(const unitig& u, vector<list_unitig>& index){
    // cout<<"add_index"<<endl;
    uint64_t size(u.sequence.size());
    if(index.size()<u.sequence.size()){
        index.resize(u.sequence.size());
    }
    // cout<<index.size()<<endl;
    index[u.sequence.size()-1].list.push_back(u);
}



uint64_t insert_fasta_in_index(const string& filename,CountMinSketch& cms,vector<list_unitig>& index, const uint64_t min_count){
    zstr::ifstream f(filename);
    string ref, useless;
    uint64_t seq_number(0);
    while (not f.eof()) {
        getline(f, useless);   // read a comment, useless
        getline(f, ref);		// read the ACGT sequence
        unitig u={ref,{filename}};
        if(not ref.empty()){
            if(cms_check_min(&cms,ref.c_str())>=min_count){
                add_index(u,index);
                seq_number++;
            }
        }
    }
    return 0;
}


struct less_than_unitig
{
    inline bool operator() (const unitig& struct1, const unitig& struct2)
    {
        return (struct1.sequence < struct2.sequence);
    }
};

uint64_t output_index(vector<list_unitig>& index, const uint64_t min_count, const string& output_filename){
    ofstream f(output_filename);
    for(uint64_t i=0;i<index.size();i++){
        if(not index[i].list.empty()){
            sort(index[i].list.begin(),index[i].list.end(),less_than_unitig());
            string output_unitig(index[i].list[0].sequence);
            vector<string> output_colors(index[i].list[0].colors);
            for(uint64_t j=1;j<index[i].list.size();j++){
                if(index[i].list[j].sequence==output_unitig){
                    output_colors.push_back(index[i].list[j].colors[0]);
                }else{
                    if(output_colors.size()>=min_count){
                        f<<">";
                        for (size_t k = 0; k < output_colors.size(); k++){
                            f<<output_colors[k]<<" ";
                        }
                        f<<"\n"<<output_unitig<<"\n";
                    }
                    output_unitig=index[i].list[j].sequence;
                    output_colors=index[i].list[j].colors;
                }
            }
            if(output_colors.size()>=min_count){
                f<<">";
                for (size_t k = 0; k < output_colors.size(); k++){
                    f<<output_colors[k]<<" ";
                }
                f<<"\n"<<output_unitig<<"\n";
            }
        }
    }
    f.close();
    return 0;
}


void help_display(){
    cout<<",.-~*´¨¯¨`*·~-.¸-This is the help message-,.-~*´¨¯¨`*·~-.,"<<endl;
    cout<<endl;
    cout<<"Main options "<<endl;
    cout<<"-p STRING CMA will process all files from this path Default: . "<<endl;
    cout<<"-c INT CMA will outpout lines seen in INT different files (or more) Default: 10"<<endl;
    cout<<"This parameter can be critical, a low value can lead to high memory usage(!)"<<endl;
    cout<<"-o STRING CMA will output shared sequences in this file Default: shared_sequences.fa "<<endl;
    cout<<"-t INT This fix the number of thread used (Not yet implemented)  Default: 1 "<<endl;
    cout<<endl;
    cout<<"Technical options (you should know what you are doing)"<<endl;
    cout<<"-w INT This fix the width of the min count sketch  Default: 1,000,000,000 "<<endl;
    cout<<"-d INT This fix the depth of the min count sketch  Default: 2 "<<endl;
}




int main(int argc, char** argv) {
    uint64_t min_count(10);
    uint64_t cms_width(1000000000);
    uint64_t cms_depth(2);
    uint64_t num_thread(1);
    string path = ".";
    string output="shared_sequences.fa";
    bool help = false;
    if (argc <2){
        help_display();
        return 0;
    }
    char c;
	while ((c = getopt (argc, argv, "o:c:w:d:p:t:h")) != -1){
		switch(c){
			case 'c':
				min_count=stoi(optarg);
				break;
            case 'w':
				cms_width=stoi(optarg);
				break;
            case 'd':
				cms_depth=stoi(optarg);
				break;
            case 't':
				num_thread=stoi(optarg);
				break;
            case 'p':
				path=(optarg);
				break;
            case 'o':
				output=(optarg);
				break;
            case 'h':
				help=true;
				break;
        }
    }
    if(help){
        help_display();
        return 0;
    }
    CountMinSketch cms;
    cms_init(&cms, (cms_width), cms_depth);
    vector<list_unitig> index;
    for (const auto & entry : std::filesystem::directory_iterator(path)){
        if(is_fasta(entry.path())){
            // cout<<"insert_fasta_in_filter "<<entry.path()<<endl;
            insert_fasta_in_filter(entry.path(),cms);
        }
    }
    for (const auto & entry : std::filesystem::directory_iterator(path)){
        if(is_fasta(entry.path())){
            // cout<<"insert_fasta_in_index "<<entry.path()<<endl;
            insert_fasta_in_index(entry.path(),cms,index,min_count);
        }
    }
    output_index(index,min_count,output);
        
    return 0;
}
