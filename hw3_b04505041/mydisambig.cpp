#include "Ngram.h"
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

double Bigram_Prob(Vocab& voc, Ngram& lm, string w1, string w2);
string Viterbi(string fw, string lw, vector<string> ZY_seq, const vector< vector<string> >& ZY_B5, const vector<string>& ZY_symbol, Vocab& voc, Ngram& lm);

int main(int argc, char *argv[]){
    if (argc != 6){	
		cerr << "Usage: ./mydisamig $input $map $lm $order $output\n";
    }
	int order = atoi(argv[4]);
	if (order != 2){
		cerr << "Only support order = 2.\n";
    }
	Vocab voc;
	Ngram lm( voc, order );
	File lmFile(argv[3], "r");
	lm.read(lmFile);
	lmFile.close();

	vector< vector<string> > ZY_B5;
	vector<string> ZYsymbol;
	ifstream ZYB5_f;
	ZYB5_f.open(argv[2]);
	if (!ZYB5_f){
		cerr << "Map file not found!\n";
    }
	string ZYB5_str;
	while (getline(ZYB5_f, ZYB5_str)) {
		if ((int)ZYB5_str.at(0) == -93) {
			ZYsymbol.push_back(ZYB5_str.substr(0,2));
			vector<string> temp_str_vect;
			ZY_B5.push_back(temp_str_vect);
		} else {
			ZY_B5.back().push_back(ZYB5_str.substr(0,2));
		}
	}
	ZYB5_f.close();

	ifstream input;
	input.open(argv[1]);
	if (!input){
		cerr << "Input text file not found!\n";
    }
	string input_str;
	ofstream output;
	output.open(argv[5]);
	int line = 1;
	while (getline(input, input_str)) {
		//cout << "Processing line " << line << endl;
		int char_count = 0;
		string output_str("<s>");
		string head_word("<s>");
		string tail_word("</s>");
		vector<string> ZY_seq;
		bool viterbi = false;
		while (char_count < input_str.size()) {			
			tail_word = head_word;
			if (input_str.at(char_count) == 32){
				char_count ++;
            }
			else {
				if(input_str.at(char_count) == -93){
						viterbi = true;
						ZY_seq.push_back(input_str.substr(char_count,2));
				}else{
					if (!viterbi){
						head_word = input_str.substr(char_count,2);	
						output_str.append(" "+head_word);
					}else {
						tail_word = input_str.substr(char_count,2);
						output_str.append(Viterbi(head_word, tail_word, ZY_seq, ZY_B5, ZYsymbol, voc, lm));
						ZY_seq.clear();
						head_word = input_str.substr(char_count,2);
						viterbi = false;
					}
				}				
				char_count += 2;
			}
		}		
		if (viterbi)
			output_str.append(Viterbi(head_word, "</s>", ZY_seq, ZY_B5, ZYsymbol, voc, lm));	
		else
			output_str.append(" </s>");
		output << output_str << "\n";
		line++;
	}
	input.close();
	output.close();
    return 0;

}

double Bigram_Prob(Vocab& voc, Ngram& lm, string w1, string w2){
	VocabIndex wid1 = voc.getIndex(w1.c_str());
	VocabIndex wid2 = voc.getIndex(w2.c_str());
	if(wid1 == Vocab_None)  //OOV
		wid1 = voc.getIndex(Vocab_Unknown);
	if(wid2 == Vocab_None)  //OOV
		wid2 = voc.getIndex(Vocab_Unknown);
	VocabIndex context[] = { wid1, Vocab_None };
	return lm.wordProb( wid2, context);
}

string Viterbi(string fw, string lw, vector<string> ZY_seq, const vector< vector<string> >& ZY_B5, const vector<string>& ZY_symbol, Vocab& voc, Ngram& lm){
	vector< vector<string> > B5_char;
	vector< vector<double> > ref_value;
	vector< vector<int> > last_ref;
	for (int ZYsb_num = 0; ZYsb_num < ZY_symbol.size(); ZYsb_num++) {
		if (ZY_seq.at(0) == ZY_symbol.at(ZYsb_num)) {
			vector<int> templastcand(ZY_B5.at(ZYsb_num).size(),0);
			vector<double> tempcandvalue(ZY_B5.at(ZYsb_num).size(),-9999);
			for (int B5_num = 0; B5_num < ZY_B5.at(ZYsb_num).size(); B5_num++) {
				templastcand.at(B5_num) = 0;
				tempcandvalue.at(B5_num) = Bigram_Prob(voc, lm, fw, ZY_B5.at(ZYsb_num).at(B5_num));
			}
			ref_value.push_back(tempcandvalue);
			last_ref.push_back(templastcand);
			B5_char.push_back(ZY_B5.at(ZYsb_num));
		}
	}
	for (int ZYseq_num = 1; ZYseq_num < ZY_seq.size(); ZYseq_num++) {
		for (int ZYsb_num = 0; ZYsb_num < ZY_symbol.size(); ZYsb_num++) {
			if (ZY_seq.at(ZYseq_num) == ZY_symbol.at(ZYsb_num)) {
				vector<int> templastcand(ZY_B5.at(ZYsb_num).size(),0);
				vector<double> tempcandvalue(ZY_B5.at(ZYsb_num).size(),-9999);
				for (int B5_num = 0; B5_num < ZY_B5.at(ZYsb_num).size(); B5_num++) {
					double max_value = -9999;
					for (int B5_char_num = 0; B5_char_num < B5_char.back().size(); B5_char_num++) {
						if (Bigram_Prob(voc, lm, B5_char.back().at(B5_char_num), ZY_B5.at(ZYsb_num).at(B5_num)) + ref_value.back().at(B5_char_num) > max_value) {
							max_value = Bigram_Prob(voc, lm, B5_char.back().at(B5_char_num), ZY_B5.at(ZYsb_num).at(B5_num)) + ref_value.back().at(B5_char_num);
							templastcand.at(B5_num) = B5_char_num;
							tempcandvalue.at(B5_num) = max_value;
						}
					}
				}
				ref_value.push_back(tempcandvalue);
				last_ref.push_back(templastcand);
				B5_char.push_back(ZY_B5.at(ZYsb_num));
			}
		}
	}
	int templastcand = 0;
	double max_value = -9999;
	for (int B5_char_num = 0; B5_char_num < B5_char.back().size(); B5_char_num++) {
		if (Bigram_Prob(voc, lm, B5_char.back().at(B5_char_num), lw) + ref_value.back().at(B5_char_num) > max_value) {
			max_value = Bigram_Prob(voc, lm, B5_char.back().at(B5_char_num), lw) + ref_value.back().at(B5_char_num);
			templastcand = B5_char_num;
		}
	}
	vector<string> tempopstr(ZY_seq.size()+1,"");
	tempopstr.back() = lw;
	
	for (int ZYseq_num = ZY_seq.size()-1; ZYseq_num >= 0; ZYseq_num--) {
		tempopstr.at(ZYseq_num) = B5_char.at(ZYseq_num).at(templastcand);
		templastcand = last_ref.at(ZYseq_num).at(templastcand);
	}
	string output_str;
	for (int B5num = 0; B5num < tempopstr.size(); B5num++)
		output_str.append(" "+tempopstr.at(B5num));
	return output_str;
}