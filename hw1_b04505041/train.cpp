#include "hmm.h"
#include <math.h>
#include <iostream>
using namespace std;

double PI[MAX_STATE];
double A[MAX_STATE][MAX_STATE];
double B[MAX_OBSERV][MAX_STATE];
double alpha[MAX_LINE][MAX_STATE];
double beta[MAX_LINE][MAX_STATE];
double gam[MAX_LINE][MAX_STATE];
double epsilon[MAX_LINE][MAX_STATE*MAX_STATE];
//void cal_alpha(int *o, int T);
//void cal_beta(int *o, int T);
//void cal_gamma(int *o, int T);
//void cal_epsilon(int *o, int T);

int main(int argc, char *argv[]){
    if (argc != 5)
        perror("Usage must be: ./train #iteration model_init.txt seq_model_01~05.txt model_01~05.txt");

    // Iteration Number
    int iternum = strtol(argv[1], NULL, 10);

    // HMM model Initialization
    HMM hmm;
    loadHMM(&hmm, argv[2]);
    const int N = hmm.state_num;
    const int O_num = hmm.observ_num;


    // Reading sequence model
    FILE *seq_model = open_or_die(argv[3], "r");
    char sample[MAX_LINE] = "";
    int sample_num = 0;
    while (fscanf(seq_model, "%s", sample) > 0)
        sample_num++;
    //printf("Number of samples: %d\n", sample_num);

    for (int i = 0; i < MAX_STATE; i++) {
        PI[i] = hmm.initial[i];
        for (int j = 0; j < MAX_STATE; j++){
            A[i][j] = hmm.transition[i][j];
        }
        for (int k = 0; k < MAX_OBSERV; k++){
            B[k][i] = hmm.observation[k][i];
        }
    }

    double gamma_all_sample[sample_num][N * (3 +O_num)];
	double epsilon_all_sample[sample_num][N * N];

    //Training
    for (int i = 1; i <= iternum; i++){
        //printf("Iteration for %d times \n", i);
        rewind(seq_model);
        int num_s = 0;
        
        // Get needed data for every observed sequences
        while (fscanf(seq_model, "%s", sample) > 0){
            int T = strlen(sample);
            int O[T];
            for (int t = 0; t < T; t++){ //A->0, B->1, C->2, D->3, E->4, F->5
                switch (sample[t]){
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

            //compute alpha
            for (int i = 0; i < N; i++){
                alpha[0][i] = PI[i] * B[O[0]][i]; // alpha_1(i) = pi_i * b_i(o_1)
            }
            for (int t = 1; t < T; t++){ // alpha_(t+1)(j) = [sigma(i=1:N)( alpha_t(i)*a_ij )] * b_j(o_(t+1))
                for (int j = 0; j < N; j++){
                    double p = 0.0;
                    for (int i = 0; i < N; i++){
                        p += alpha[t - 1][i] * A[i][j];
                    }
                    alpha[t][j] = p * B[O[t]][j];
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
                        beta[t][i]+= A[i][j] * B[O[t + 1]][j] * beta[t + 1][j];
                    }
                }
            }


            double prob_sample = 0;
            double gamma_T_sum[N];    // sigma(t=1:T-1)gamma_t(i)
		    double gamma_T1_sum[N];    // sigma(t=1:T)gamma_t(i)
		    double gamma_observ_sum[N * O_num];  // sigma(t=1:T, o_t=v_k)gamma_t(i)
		    double gamma_one_sample[N * (3 + O_num)];

            for ( int i = 0; i < N; i++ ) {
				gamma_T1_sum[i] = 0;
				gamma_T_sum[i] = 0;
				prob_sample += alpha[T/2][i] * beta[T/2][i];
			}
			for ( int i = 0; i < N * O_num; i++ ){
				gamma_observ_sum[i] = 0;
            }

            //compute gamma
            for (int t = 0; t < T; t++){
                for (int i = 0; i < N; i++){
                    gam[t][i] = alpha[t][i] * beta[t][i] / prob_sample; //gamma_t(i) = alpha_t(i)*beta_t(i) / sigma(i=1:6)( alpha_t(i)*beta_t(i) )
                    gamma_observ_sum[O[t]*N + i] += gam[t][i];      // sigma(t=1:T, o_t=v_k)gamma_t(i)
                    gamma_T_sum[i] += gam[t][i];    //sigma(t=1:T)gamma_t(i)
                    if(t != T-1){   
                        gamma_T1_sum[i] += gam[t][i]; // sigma(t=1:T-1)gamma_t(i)
                    }
                }    
            }
            for ( int i = 0; i < N; i++ ) {
				gamma_one_sample[i] = gam[0][i];
				gamma_one_sample[N + i] = gamma_T1_sum[i];
				gamma_one_sample[2 * N + i] = gamma_T_sum[i];
			}
            for ( int i = 0; i < N * O_num; i++ )
				gamma_one_sample[3*N + i] = gamma_observ_sum[i];

            //compute epsilon
            double epsilon_one_sample[N*N];
            for(int i = 0 ; i < N*N ; i++){
                epsilon_one_sample[i] = 0;
            }
            
            // epsilon_t(i,j) = alpha_t(i)*a_ij*b_j(o_(t+1))*beta_(t+1)(j) / sigma(i=1:6)sigma(j=1:6)( alpha_t(i)*a_ij*b_j(o_(t+1))*beta_(t+1)(j) )
            for (int t = 0; t < T - 1; t++){
                for (int i = 0; i < N; i++){
                    for (int j = 0; j < N; j++){
                        epsilon[t][i * N + j]  = alpha[t][i] * A[i][j] * B[O[t + 1]][j] * beta[t + 1][j] / prob_sample;
                        epsilon_one_sample[i * N + j] += epsilon[t][i * N + j]; // sigma(t=1:T-1)epsilon_t(i,j)
                    }
                }
            }

            for ( int i = 0; i < N * (3 + O_num); i++ )
				gamma_all_sample[num_s][i] = gamma_one_sample[i];
			for ( int i = 0; i < N * N; i++ ){
				epsilon_all_sample[num_s][i] = epsilon_one_sample[i];
            }

            num_s++;


        }

        //Update new data
        double new_A[N*N];
		double new_B[N*O_num];
        double sum_gamma_T1_sum[N];
		double sum_gamma_T_sum[N];

        for (int i = 0; i < N; i++){ //update new PI
            PI[i] = 0;
            for(int s = 0 ; s < sample_num ; s++){
                PI[i] += gamma_all_sample[s][i];
            }
            PI[i] /= sample_num;
        }
        for (int i = 0; i < N; i++){ // update new A
            sum_gamma_T1_sum[i] = 0;
            for(int s = 0 ; s < sample_num ; s++){
                sum_gamma_T1_sum[i] += gamma_all_sample[s][N + i];  // sum of sigma(t=1:T-1)gamma_t(i): sigma sigma(t=1:T-1)gamma_t(i) 
            }
            for (int j = 0 ; j < N; j++){
                new_A[i * N + j] = 0;
                for(int s = 0 ; s < sample_num ; s++){
                    new_A[i * N + j] +=  epsilon_all_sample[s][i * N + j];
                }
                new_A[i * N + j] /= sum_gamma_T1_sum[i];
            }
        }

        for (int i = 0; i < N; i++){ //update new B
            sum_gamma_T_sum[i] = 0;
            for(int s = 0 ; s < sample_num ; s++){
                sum_gamma_T_sum[i] += gamma_all_sample[s][2 * N + i];  // sum of sigma(t=1:T)gamma_t(i): sigma sigma(t=1:T)gamma_t(i) 
            }
            for (int k = 0; k < 6; k++){
                new_B[k * N + i] = 0;
                for(int s = 0 ; s < sample_num ; s++){
                    new_B[k * N + i] += gamma_all_sample[s][3*N + k*O_num + i];
                }
                new_B[k * N + i] /= sum_gamma_T_sum[i];
            }
        }

        for( int i = 0 ; i < N; i++ ) {
		    hmm.initial[i] = PI[i];
			for( int j = 0 ; j < N; j++ ){
                A[i][j] = new_A[i * N + j];
		        hmm.transition[i][j] = A[i][j];
            }
			for( int k = 0 ; k < O_num; k++ ){
                B[k][i] = new_B[k * O_num +i];
				hmm.observation[k][i] = B[k][i];
            }
		}

        //show updated model every 100 iterations
        if (i % 100 == 0)
        {
            printf("\nAfter iterations %d\n\n", i);
            printf("\ninitial: %d\n", N);
            for (int i = 0; i < N; i++)
                printf("%.5lf ", hmm.initial[i]);

            printf("\n\ntransition: %d\n", N);
            for (int i = 0; i < N; i++)
            {
                for (int j = 0; j < N; j++)
                    printf("%.5lf ", hmm.transition[i][j]);
                printf("\n");
            }
            printf("\nobservation: %d\n", O_num);
            for (int k = 0; k < O_num; k++)
            {
                for (int j = 0; j < N; j++)
                    printf("%.5lf ", hmm.observation[k][j]);
                printf("\n");
            }
            printf("\n\n");
        }
    }
    FILE *final_model = open_or_die(argv[4], "w");
    // Save the trained model
    dumpHMM(final_model, &hmm);
    return 0;
}

// void cal_alpha(int *o, int T)
// {
//     for (int i = 0; i < N; i++)
//     {
//         alpha[0][i] = PI[i] * B[o[0]][i]; // alpha_1(i) = pi_i * b_i(o_1)
//     }
//     for (int t = 1; t < T; t++)
//     { // alpha_(t+1)(j) = [sigma(i=1:N)( alpha_t(i)*a_ij )] * b_j(o_(t+1))
//         for (int j = 0; j < N; j++)
//         {
//             double p = 0.0;
//             for (int i = 0; i < N; i++)
//             {
//                 p += alpha[t - 1][i] * A[i][j];
//             }
//             alpha[t][j] = p * B[o[t]][j];
//         }
//     }
// }
// void cal_beta(int *o, int T)
// {
//     for (int i = 0; i < N; i++)
//     { // beta_T(i) = 1
//         beta[T][i] = 1;
//     }
//     for (int t = t - 2; t >= 0; t--)
//     { // beta_(t)(i) = sigma(j=1:6)( a_ij*b_j(o_(t+1))*beta_(t+1)(j))
//         for (int i = 0; i < N; i++)
//         {
//             double p = 0.0;
//             for (int j = 0; j < N; j++)
//             {
//                 p += A[i][j] * B[o[t + 1]][j] * beta[t + 1][j];
//             }
//             beta[t][i] = p;
//         }
//     }
// }
// void cal_gamma(int *o, int T)
// {
//     double temp_gamma[MAX_STATE];
//     for (int t = 0; t < T; t++)
//     {
//         double temp = 0.0;
//         for (int i = 0; i < N; i++)
//         {
//             temp_gamma[i] = alpha[t][i] * beta[t][i];
//             temp += temp_gamma[i];
//         }
//         for (int i = 0; i < N; i++)
//         {
//             temp_gamma[i] /= temp;
//             gam[i] += temp_gamma[i];
//             gam_full[i] += temp_gamma[i];
//             gam_bucket[o[t]][i] += temp_gamma[i];
//             if (t == 0)
//                 temp_pi[i] += temp_gamma[i];
//         }
//     }
//     for (int i = 0; i < N; i++)
//     {
//         gam[i] -= temp_gamma[i];
//     }
//     return;
// }
// void cal_epsilon(int *o, int T)
// {
//     double temp_epsilon[MAX_STATE][MAX_STATE];
//     // epsilon_t(i,j) = alpha_t(i)*a_ij*b_j(o_(t+1))*beta_(t+1)(j) / sigma(i=1:6)sigma(j=1:6)( alpha_t(i)*a_ij*b_j(o_(t+1))*beta_(t+1)(j) )
//     for (int t = 0; t < T - 1; t++)
//     {
//         double temp = 0.0;
//         for (int i = 0; i < N; i++)
//         {
//             for (int j = 0; j < N; j++)
//             {
//                 temp_epsilon[i][i] = alpha[t][i] * A[i][j] * B[o[t + 1]][j] * beta[t + 1][j];
//                 temp += temp_epsilon[i][j]; //sigma(i=1:6)sigma(j=1:6)( alpha_t(i)*a_ij*b_j(o_(t+1))*beta_(t+1)(j) )
//             }
//         }
//         for (int i = 0; i < N; i++)
//         {
//             for (int j = 0; j < N; j++)
//             {
//                 epsilon[i][j] += (temp_epsilon[i][j] / temp); // summing epsilon
//             }
//         }
//     }
// }
