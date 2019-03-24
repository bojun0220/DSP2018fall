#include "hmm.h"
#include <math.h>
#include <iostream>
#include <cstring>
#include <string>

using namespace std;

double PI[5][MAX_STATE];
double A[5][MAX_STATE][MAX_STATE];
double B[5][MAX_OBSERV][MAX_STATE];
double alpha[MAX_LINE][MAX_STATE];
double beta[MAX_LINE][MAX_STATE];
double delta[MAX_LINE][MAX_STATE];
int main(int argc , char *argv[]){
    if (argc != 4) 
		perror ("Usage must be: ./test modellist.txt testing_data.txt result.txt");

    // Loading Trained Model
	HMM hmm[5];
	load_models( argv[1], hmm, 5);

	// Get model list	
	FILE *modellist = open_or_die( argv[1], "r");
	char modelname[MAX_LINE] = "";
	char modelname_list[5][MAX_LINE];
	int model_number = 0;
	while ( fscanf( modellist, "%s", modelname) > 0 ) {
		strcpy(modelname_list[model_number], modelname);
		model_number++;
	}

	// Loading testing sequences, getting the number of sequences
	FILE *testing_data = open_or_die( argv[2], "r");
	char test_seq[MAX_LINE] = "";
	int seq_num = 0;
	const int N = hmm[0].state_num;
	const int O_num = hmm[0].observ_num;
	while ( fscanf( testing_data, "%s", test_seq) > 0 ){
		seq_num++;
	}
	//printf("Number of sequences: %d\n", seq_num);

	for(int m = 0 ; m < 5 ; m++){	//Lodaing 5 models of initail , transition and observation into PI,A,B
		for (int i = 0; i < MAX_STATE; i++) {
        	PI[m][i] = hmm[m].initial[i];
        	for (int j = 0; j < MAX_STATE; j++){
           	 	A[m][i][j] = hmm[m].transition[i][j];
        	}
        	for (int k = 0; k < MAX_OBSERV; k++){
            	B[m][k][i] = hmm[m].observation[k][i];
        	}	
    	}
	}

	rewind(testing_data);
	char** testing_answer = (char**)malloc( sizeof(char*) * seq_num );   
	double prob_corresponding[seq_num];  
	int ns = 0;
	model_number = 0;

	while ( fscanf( testing_data, "%s", test_seq) > 0 ){
		int T = strlen(test_seq);
		testing_answer[ns] = (char*)malloc( sizeof(char) * MAX_LINE );
		int O[T];
		for (int t = 0; t < T; t++){ //A->0, B->1, C->2, D->3, E->4, F->5
			switch (test_seq[t]){
                case 'A' :
                    O[t] = 0;
                    break;
                case 'B' :
                    O[t] = 1;
                        break;
                case 'C' :
                    O[t] = 2;
                        break;
                case 'D' :
                    O[t] = 3;
                        break;
                case 'E' :
          		    O[t] = 4;
                        break;
                case 'F' :
                    O[t] = 5;
                        break;
            }
        }
		double prob_star[5];
		double max = 0;

		
		for(int m = 0 ; m < 5 ; m ++){
			//compute alpha
			for (int i = 0; i < N; i++){
                alpha[0][i] = PI[m][i] * B[m][O[0]][i]; // alpha_1(i) = pi_i * b_i(o_1)
            }
            for (int t = 1; t < T; t++){ // alpha_(t+1)(j) = [sigma(i=1:N)( alpha_t(i)*a_ij )] * b_j(o_(t+1))
                for (int j = 0; j < N; j++){
                    double p = 0.0;
                    for (int i = 0; i < N; i++){
                        p += alpha[t - 1][i] * A[m][i][j];
                    }
                    alpha[t][j] = p * B[m][O[t]][j];
                }
            }

			//compute beta
			for (int i = 0; i < N; i++){ // beta_T(i) = 1
                beta[T-1][i] = 1;
            }
            for (int t = T - 2; t >= 0; t--){ // beta_(t)(i) = sigma(j=1:6)( a_ij*b_j(o_(t+1))*beta_(t+1)(j))
                for (int i = 0; i < N; i++){
                    beta[t][i]= 0.0;
                    for (int j = 0; j < N; j++){
                        beta[t][i]+= A[m][i][j] * B[m][O[t + 1]][j] * beta[t + 1][j];
                    }
                }
            }

			//compute delta
			for(int i = 0 ; i < N ; i++){	// delta_1(i) = pi_i * b_i(o_1)
				delta[0][i] = PI[m][i] * B[m][O[0]][i];
			}
			for(int t = 1 ; t < T ; t++){	//delta_t(j) = max(i=1:N)( delta_(t-1)(i)*a_ij )*b_j(o_t)
				for(int j = 0 ; j < N ; j++){
					max = 0;
					for(int i = 0 ; i < N ; i++){
						if (max < delta[t-1][i] * A[m][i][j]){
							max = delta[t-1][i] * A[m][i][j]; // max(i=1:N) delta(t-1)i*aij
						}
					}
					delta[t][j] = max * B[m][O[t]][j];
				}
			}

			prob_star[m] = 0;
			for(int i = 0 ; i < N ;i++){
				if(prob_star[m] < delta[T-1][i]){
					prob_star[m] = delta[T-1][i];
				}

			}
		}
		prob_corresponding[ns] = 0;
		for (int m = 0; m < 5; m++ ) {
			if ( prob_corresponding[ns] < prob_star[m] ) {
				prob_corresponding[ns] = prob_star[m];
				model_number = m;
			}
		}
		testing_answer[ns] = modelname_list[model_number];
		ns++;
	}
	FILE *test_result = open_or_die( argv[3], "w");
	ns = 0;
	for ( ns = 0; ns < seq_num; ns++ ){
		fprintf( test_result, "%s %.6e\n", testing_answer[ns], prob_corresponding[ns]);
	}
	//Just to compare the answer with the given "result data.txt",not necessary
	if ( strcmp(argv[2], "testing_data1.txt") == 0 ) {	
		ns = 0;
		FILE *test_answer = open_or_die( "testing_answer.txt", "r");
		char** true_answer = (char**)malloc( sizeof(char*) * seq_num );
		for (int n = 0; n < seq_num; n++ ){
			true_answer[n] =  (char*)malloc( sizeof(char) * MAX_LINE );
		}
		while ( fscanf( test_answer, "%s", modelname) > 0 ) {
			strcpy(true_answer[ns], modelname);
			ns++;
		}	
		int num_correct = 0;
		for (int n = 0; n < seq_num; n++ ) {
			if ( strcmp(true_answer[n], testing_answer[n]) == 0 )
				num_correct++;
		}
		double accuracy = (double)num_correct/(double)seq_num; 
		printf("accurracy: %.5lf\n", accuracy);
		FILE *acc = open_or_die( "acc.txt", "w");
		fprintf( acc, "%.5lf", accuracy);
	}
	return 0;





}