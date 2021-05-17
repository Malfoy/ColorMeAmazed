#include <vector>
#include <string>
#include <filesystem>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <unistd.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "strict_fstream.hpp"
#include "zstr.hpp"
#include "count_min_sketch.h"



//TODO multithread



using namespace std;



struct unitig{
    string sequence;
    vector<uint32_t>colors;
};



struct list_unitig{
    vector<unitig> list;
};



struct index_color{
    vector<list_unitig> size_list;
    vector<string> filenames;
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
            seq_number++;
            cms_add(&cms,ref.c_str());
            // cout<<"i add:   "<<ref<<endl;
        }
        ref=useless="";
        
    }
    return seq_number;
}



void add_index(const unitig& u, vector<list_unitig>& index){
    // cout<<"add_index"<<endl;
    if(index.size()<u.sequence.size()){
        index.resize(u.sequence.size());
    }
    // cout<<index.size()<<endl;
    index[u.sequence.size()-1].list.push_back(u);
}



uint64_t insert_fasta_in_index(const string& filename,CountMinSketch& cms, index_color& index, const uint64_t min_count){
    zstr::ifstream f(filename);
    index.filenames.push_back(filename);
    uint32_t filename_indice(index.filenames.size()-1);
    string ref, useless;
    uint64_t seq_number(0);
    while (not f.eof()) {
        getline(f, useless);   // read a comment, useless
        getline(f, ref);		// read the ACGT sequence
        unitig u={ref,{filename_indice}};
        if(not ref.empty()){
            if((uint64_t)cms_check_min(&cms,ref.c_str())>=min_count){
                add_index(u,index.size_list);
                seq_number++;
            }
            ref.clear();
        }
    }
    return seq_number;
}


struct less_than_unitig
{
    inline bool operator() (const unitig& struct1, const unitig& struct2)
    {
        return (struct1.sequence < struct2.sequence);
    }
};



uint64_t output_index(index_color& index, const uint64_t min_count, const string& output_filename,uint64_t& overestimated_sequence){
    uint64_t output_sequence(0);
    ofstream f(output_filename);
    for(uint64_t i=0;i<index.size_list.size();i++){
        if(not index.size_list[i].list.empty()){
            sort(index.size_list[i].list.begin(),index.size_list[i].list.end(),less_than_unitig());
            string output_unitig(index.size_list[i].list[0].sequence);
            vector<uint32_t> output_colors(index.size_list[i].list[0].colors);
            for(uint64_t j=1;j<index.size_list[i].list.size();j++){
                if(index.size_list[i].list[j].sequence==output_unitig){
                    output_colors.push_back(index.size_list[i].list[j].colors[0]);
                }else{
                    if(output_colors.size()>=min_count){
                        f<<">";
                        for (size_t k = 0; k < output_colors.size(); k++){
                            f<<index.filenames[output_colors[k]]<<" ";
                        }
                        f<<"\n"<<output_unitig<<"\n";
                        output_sequence++;
                    }else{
                        overestimated_sequence++;
                    }
                    output_unitig=index.size_list[i].list[j].sequence;
                    output_colors=index.size_list[i].list[j].colors;
                }
            }
            if(output_colors.size()>=min_count){
                f<<">";
                for (size_t k = 0; k < output_colors.size(); k++){
                    f<<index.filenames[output_colors[k]]<<" ";
                }
                f<<"\n"<<output_unitig<<"\n";
                output_sequence++;
            }
        }
    }
    f.close();
    return output_sequence;
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


uint64_t getMemorySelfMaxUsed (){
	uint64_t result = 0;
	struct rusage usage;
	if (getrusage(RUSAGE_SELF, &usage)==0)  {  result = usage.ru_maxrss;  }
	return result/1024;
}

string intToString(uint64_t n){
	string result;
	uint64_t order=1000000000000000000;
	bool started(false);
	while(order!=1){
		if(n/order>=1){
			string local( to_string(n/order));
			if(started){
				if(local.size()==2){
					result+='0';
				}
				if(local.size()==1){
					result+="00";
				}
			}
			result+=local+",";
			started=true;
			n%=order;
		}else if (started){
			result+="000,";
		}
		order/=1000;
	}
	string local( to_string(n));
	if(started){
		if(local.size()==2){
			result+='0';
		}
		if(local.size()==1){
			result+="00";
		}
	}
	result+=local;

	return result;
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
				num_thread+=stoi(optarg);
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
    uint64_t file_read(0);
    uint64_t sequence_read(0);
    uint64_t sequence_indexed(0);
    uint64_t output_sequence(0);
    uint64_t overestimated_sequence(0);

    CountMinSketch cms;
    cms_init(&cms, (cms_width), cms_depth);
    cout<<"Insertion in the filter"<<endl;
    auto start = std::chrono::system_clock::now();

    for (const auto & entry : std::filesystem::directory_iterator(path)){
        if(is_fasta(entry.path())){
            sequence_read+=insert_fasta_in_filter(entry.path(),cms);
            file_read++;
            if(file_read%10==0)
                cout<<"-"<<flush;
        }
    }
    cout<<endl;
    auto middle = std::chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = middle - start;
	cout << "Elapsed time: " << elapsed_seconds.count() << "s ("<< (double)elapsed_seconds.count()/file_read<<"s per file)"<<endl;
    cout<<"Memory used: " <<intToString(getMemorySelfMaxUsed())<<" MB ("<<intToString(getMemorySelfMaxUsed()*1024*1024/file_read)<<" B per file)"<<endl;

    cout<<"Frequent sequences indexing"<<endl;
    index_color index;
    file_read=0;
    for (const auto & entry : std::filesystem::directory_iterator(path)){
        if(is_fasta(entry.path())){
            sequence_indexed+=insert_fasta_in_index(entry.path(),cms,index,min_count);
            file_read++;
            if(file_read%10==0)
                cout<<"-"<<flush;
        }
    }
    cout<<endl;
     auto end = std::chrono::system_clock::now();
    elapsed_seconds = end- middle ;
	cout << "Elapsed time: " << elapsed_seconds.count() << "s ("<< (double)elapsed_seconds.count()/file_read<<"s per file)"<<endl;
    cout<<"Memory used: " <<intToString(getMemorySelfMaxUsed())<<" MB ("<<intToString(getMemorySelfMaxUsed()*1024*1024/file_read)<<" B per file)"<<endl;
    cout<<"Output shared sequences"<<endl;
    output_sequence=output_index(index,min_count,output,overestimated_sequence);

    cout<<"I read "<<intToString(file_read)<<" files containing "<<intToString(sequence_read)<<" sequences" <<endl;
    cout<<intToString(sequence_indexed)<<" sequences passed the filter and indexed"<<endl;
    cout<<intToString(overestimated_sequence)<<" sequences were overestimated"<<endl;
    cout<<intToString(output_sequence)<<" sequences were actually output"<<endl;
    
        
    return 0;
}
