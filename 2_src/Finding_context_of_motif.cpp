/*#########################################
## 04/05/2023
## Par Elyna Bouchereau
## Fichier: Finding_context_of_motif.cpp
###########################################*/
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

struct sequence {
	string header;
	string sequence;
};

void fasta_file_print(string file_name){
	sequence seq;
	ifstream fasta_file(file_name);
	string line;
	
	while(getline(fasta_file, line)){
		if(line.empty()) continue;
		if(line[0] == '>'){
			if(!seq.header.empty()) cout << seq.header << " : " << seq.sequence << endl;
			seq.header = line.substr(1);
			seq.sequence.clear();
		} else {
			seq.sequence += line;
		}
	}
	if(!seq.header.empty()) cout << seq.header << " : " << seq.sequence << endl;

}

void fasta_handler(string file_name, vector<sequence> &seq_list){
	sequence seq;
	ifstream fasta_file(file_name);
	string line;

	//cout << "Fasta handler: ";
	while(getline(fasta_file, line)){
		if(line.empty()) continue;
		if(line[0] == '>'){
			if(!seq.header.empty()) seq_list.push_back(seq);
			seq.header = line.substr(1);
			seq.sequence.clear();
		} else {
			seq.sequence += line;
		}
	}
	if(!seq.header.empty()) seq_list.push_back(seq);
	//cout << "DONE\n";
}

void fasta_print(vector<sequence> &seq_list){
	for(int i=0;i<seq_list.size();i++){
		cout << seq_list[i].header << " : " << seq_list[i].sequence << endl;
	}
}

void find_motif(vector<sequence> &seq_list, string motif){
	size_t position;

	for(int i=0;i<seq_list.size();i++){
		position = seq_list[i].sequence.find(motif);
		if (position != string::npos) {
			cout << position + 20 << " " << seq_list[i].header.substr(0,10) << endl;
		}
	}
}


//################################################################################## Doh, Calligraphy, Broadway, Fraktur, Poison
/*                                                                                                                   
     *****   **    **           **                *****  *      ***** *     **    
  ******  ***** *****        *****             ******  *     ******  **    **** * 
 **   *  *  ***** *****     *  ***            **   *  *     **   *  * **    ****  
*    *  *   * **  * **         ***           *    *  *     *    *  *  **    * *   
    *  *    *     *           *  **              *  *          *  *    **   *     
   ** **    *     *           *  **             ** **         ** **    **   *     
   ** **    *     *          *    **            ** **         ** **     **  *     
   ** **    *     *          *    **          **** **         ** **     **  *     
   ** **    *     *         *      **        * *** **         ** **      ** *     
   ** **    *     **        *********           ** **         ** **      ** *     
   *  **    *     **       *        **     **   ** **         *  **       ***     
      *     *      **      *        **    ***   *  *             *        ***     
  ****      *      **     *****      **    ***    *          ****          **     
 *  *****           **   *   ****    ** *   ******          *  *****              
*     **                *     **      **      ***          *     **               
*                       *                                  *                      
 **                      **                                 **                    
*/
int main(int argc, char *argv[]){
	//*####			DEFINE VARIABLES		#######
	string genome_file;
	vector<string> context_seq;
	vector<sequence> seq_list(0);

	//*####			CHARGE PARAMETERS		#######
	genome_file = argv[1];
	for(int i=1; i < argc; i++){
		context_seq.push_back(argv[i]);
	}

	//*####			EXECUTION MAIN			#######
	fasta_handler(genome_file, seq_list);
	//fasta_print(seq_list);
	for(int j=0;j<context_seq.size();j++){
		find_motif(seq_list, context_seq[j]);
	}
}