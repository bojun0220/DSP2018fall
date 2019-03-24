#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
using namespace std;
#define ZY_num 37
vector<string> split(string str_to_split, char delimeter){//string split implementation by using delimiter as a character.
    stringstream ss(str_to_split);
    string item;
    vector<string> splittedString;
    while (getline(ss, item, delimeter)){
       splittedString.push_back(item);
    }
    return splittedString;
}
int main(){
    ifstream B5ZY_f;//Read Big5-ZhuYin.map
	B5ZY_f.open("Big5-ZhuYin.map");
	if (!B5ZY_f){
		cerr << "The file 'Big5-ZhuYin.map' is not found!\n";
	}
	string B5ZY_str;
	vector<string> ZY[ZY_num] ;
	for(int i = 0 ; i < ZY_num ; i++){
		ZY[i].push_back(" ");
	}
	while (getline(B5ZY_f, B5ZY_str)) {
		vector<string> data = split(B5ZY_str , ' ');
		string Chi = data[0]; //The Chinese Character
		vector<string> zy_list = split(data[1] , '/'); //the corresponding list of Zhu-Yin notations
		for(int i = 0 ; i < zy_list.size() ; i++){ //establish Zhu_Yin to BIG5 list
			string zy_test = zy_list[i].substr(0,2);
			for(int j = 0 ; j < ZY_num ; j++){
				if(ZY[j].front() == zy_test){
					if(find(ZY[j].begin(), ZY[j].end(), Chi) == ZY[j].end()){
						ZY[j].push_back(Chi);}
					break;
				}else if(ZY[j].front() == " "){
					ZY[j].erase(ZY[j].begin());
					ZY[j].insert(ZY[j].begin(), zy_test);
					ZY[j].push_back(Chi);
					break;
				}
			}
		}	
	}
	B5ZY_f.close();
	ofstream ZYB5_f;
	ZYB5_f.open("ZhuYin-Big5.map");
	for (int i = 0; i < ZY_num ; i++) {
		ZYB5_f << ZY[i].front() << "  ";
		for(int j = 1 ; j < ZY[i].size() ; j++){
			ZYB5_f << ZY[i].at(j) << " ";
		}
		ZYB5_f <<"\n";
		for (int k = 1 ; k < ZY[i].size() ; k++){
			ZYB5_f << ZY[i].at(k) << "  " << ZY[i].at(k) << "\n";
		}
	}
	return 0;
}