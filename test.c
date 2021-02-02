#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common.h"
#include "header.h"
#include "gamma.h"


/*Claim Functions defined in this file*/
void input1(FILE *infile);
void input2(FILE *infile);
double CSDapprx_scaling(long* TYPE, long* newseq, long nbranch, long ntype, const double RHO);

long Ne; /*effective population size*/
double NumLineage; /*the number of lineages to stop */
double THETA; /*mutation rate*/
double RHO; /*recombination rate*/
double SEQLEN; /*sequence length*/

long M; /*number of genes*/
long L; /*number of loci*/
double Mu; /*stationary distn of allele 0 in the samples*/
long *TYPE_copy; /*store observed haplotypes*/
long ntype_copy;
double P[2][2]; /*mutation matrix*/
long *new; /*new sequence */
double *positions; /*genetic position for each site */

double adj_rate;
double gamma2[2][NODE_MAX];


int main(int argc, char *argv[]){
    FILE *in1, *in2;
    double logCSD;

    fopen_s(&in1, "C:\\Users\\blabl\\Desktop\\tests\\test1","r");
    input1(in1);
    fclose(in1);

    fopen_s(&in2, "C:\\Users\\blabl\\Desktop\\tests\\test2","r");
    input2(in2);
    fclose(in2);

    //THETA = 4.0 * Ne * THETA * SEQLEN / (double)L;
    RHO = 4.0 * Ne * RHO;

    logCSD = CSDapprx_scaling(TYPE_copy, new, M, ntype_copy, RHO);

    fprintf(stderr, "logCSD = %lf \t", logCSD);

    free(TYPE_copy);
    free(positions);
    free(new);

    return(0);

}


double CSDapprx_scaling(long* TYPE, long* newseq, long nbranch, long ntype, const double RHO)
{
	long i, j, * index, index_len = 0, d, a;
	int e;
	double pj, alpha_sum, * prop, * alpha_prime, * alpha_hat, logwt_scale; 

	index = long_vec_init(L);
	

	/*set up index-loci ancestral in new*/
	for (i = 0; i < L; i++) {
		if (*(newseq + i) >= 0) {
			*(index + index_len) = i;
			index_len++;
		}
	}


	alpha_prime = dou_vec_init(ntype);


	/*calculate prop and gamma */
	/*prop is the proportion of a chromosome in the current configuration
	that is of type 0 at each site. It's a vector of length index_len */
	prop = dou_vec_init(index_len);
	calcProp(prop, TYPE, ntype, index, index_len);


    gammaMat2(gamma2, THETA);
	alpha_hat = dou_vec_init(ntype);

	/*calculate alpha_1(x), for x = 1,...,k. alpha is of length ntype */
	d = newseq[index[0]];
	if (nbranch >= NODE_MAX) a = NODE_MAX - 1;
	else a = nbranch;
	for (i = 0; i < ntype; i++) {
		if (*(TYPE + i * (L + 1) + L) > 0) {
			e = *(TYPE + i * (L + 1) + index[0]);
			if (e >= 0) {
				if (d == e) {
					alpha_prime[i] = gamma2[d][a] / (double)nbranch;
				}
				else {
					alpha_prime[i] = (1.0 - gamma2[d][a]) / (double)nbranch;
				}
			}
			else {
				if (d == 0) {
					alpha_prime[i] = (prop[0] * gamma2[d][a] + (1.0 - prop[0]) * (1.0 - gamma2[d][a])) / (double)nbranch;
				}
				else {
					alpha_prime[i] = ((1.0 - prop[0]) * gamma2[d][a] + prop[0] * (1.0 - gamma2[d][a])) / (double)nbranch;
				}
			}
			
		}
	}
	alpha_sum = 0.0;
	for (i = 0; i < ntype; i++) {
		if (*(TYPE + i * (L + 1) + L) > 0) {
			alpha_sum += *(TYPE + (L + 1) * i + L) * alpha_prime[i];
		}
	}
	for (i = 0; i < ntype; i++) {
		alpha_hat[i] = alpha_prime[i] / alpha_sum;
	}
	logwt_scale = log(alpha_sum);

	/*calculate alpha_j, for j=2,...,S */
	for (j = 1; j < index_len; j++) {
		/*calculate recombination rate, p_j, in (A5) */
		pj = exp(-1.0 * (positions[*(index + j)] - positions[index[j - 1]]) * adj_rate * RHO / (double)nbranch);
		if (nbranch >= NODE_MAX) a = NODE_MAX - 1;
		else a = nbranch;

		d = newseq[index[j]];
		if (d != 0 && d != 1) fprintf(stderr, "newseq error!\n");
		for (i = 0; i < ntype; i++) {
			if (*(TYPE + i * (L + 1) + L) > 0) {
				e = *(TYPE + i * (L + 1) + index[j]);

				if (e >= 0) {
					if (d == e) {
						alpha_prime[i] = (pj * alpha_hat[i] + (1.0 - pj) / (double)nbranch) * gamma2[d][a];
					
					}
					else {
						alpha_prime[i] = (pj * alpha_hat[i] + (1.0 - pj) / (double)nbranch) * (1.0 - gamma2[d][a]);
						
					}
				}
				else {
					if (d == 0) {
						alpha_prime[i] = (pj * alpha_hat[i] + (1.0 - pj) / (double)nbranch) * (prop[j] * gamma2[d][a] + (1.0 - prop[j]) * (1.0 - gamma2[d][a]));
						

					}
					else {
						alpha_prime[i] = (pj * alpha_hat[i] + (1.0 - pj) / (double)nbranch) * ((1.0 - prop[j]) * gamma2[d][a] + prop[j] * (1.0 - gamma2[d][a]));
						
					}
				}

			}
			//if (newalpha[i] <= 0)
				//fprintf(stderr, "Error! alpha[%ld] <= 0!\n", i);
		}
		alpha_sum = 0.0;
		for (i = 0; i < ntype; i++) {
			if (*(TYPE + i * (L + 1) + L) > 0) {
				alpha_sum += *(TYPE + (L + 1) * i + L) * alpha_prime[i];
			}
		}
		for (i = 0; i < ntype; i++) {
			alpha_hat[i] = alpha_prime[i] / alpha_sum;
		}
		logwt_scale += log(alpha_sum);
	}

	free(index);
	free(prop);
	free(alpha_prime);
	free(alpha_hat);


	return(logwt_scale);
}







void input1(FILE *infile)
{
    char string[200];
    char iden[2];
    short lines;
    short flag=1;
    short i, j;
    
	Ne = -1;
    RHO = -1; /*recombination rate*/
    THETA = -1; /*mutation rate*/
	SEQLEN = -1;
	NumLineage = -1;
    lines = 0;

    fprintf(stderr,"\n Data in file 1: \n");
    
    while(flag){
        if(fgets(string,sizeof(string),infile)==NULL){
            fprintf(stderr, "Error in file \n");
            abort();
        }
        else{
            i=0;
            while(string[i]==' ' || string[i]=='\t') i++;
            iden[0]=string[i];
            switch(iden[0]){
                case EOF:
                    flag=0;
                    break;
                case '#':
                    flag=0;
                    break;
                case '\n':
                    break;
        

				case 'E':
					if (Ne != -1) {
						fprintf(stderr, "Effective population size redefined on line %d \n", lines);
						break;
					}
					else {
						while (string[i] != '=') i++;
						i++;
						sscanf_s(string + i, "%ld", &Ne);
						fprintf(stderr, "effective population size %ld \n", Ne);
						break;
					}

                
                case'R':
                    if(RHO != -1){
                        fprintf(stderr,"Recombination rate redefined on line %d. Ignored \n",lines);
                        break;
                    }
                    else{
                        while(string[i]!='=') i++;
                        i++;
						sscanf_s(string+i,"%lf",&RHO);
                        fprintf(stderr,"recombination rate %3.9lf per generation per unit physical distance\n",RHO);
                        break;
                    }
                    
                case'M':
                    if(THETA != -1){
                        fprintf(stderr,"Mutation rate redefined on line %d. Ignored \n",lines);
                        break;
                    }
                    else{
                        while(string[i]!='=') i++;
                        i++;
						sscanf_s(string+i,"%lf",&THETA);
                        fprintf(stderr,"mutation rate %3.9lf per generation per base pair\n",THETA);
                        break;
                    }

				case'L':
					if (SEQLEN != -1) {
						fprintf(stderr, "Sequence length redefined on line %d. Ignored \n", lines);
						break;
					}
					else {
						while (string[i] != '=') i++;
						i++;
						sscanf_s(string + i, "%lf", &SEQLEN);
						fprintf(stderr, "sequence length %lf \n", SEQLEN);
						break;
					}

				case'S':
					if (NumLineage != -1) {
						fprintf(stderr, "Number of lineages redefined on line %d. Ignored \n", lines);
						break;
					}
					else {
						while (string[i] != '=') i++;
						i++;
						sscanf_s(string + i, "%lf", &NumLineage);
						fprintf(stderr, "number of lineages %lf \n", NumLineage);
						break;
					}


                default:
                    fprintf(stderr,"Non-standard input line %d in file 1. Ignored",lines);
                    break;
            }
            lines++;
        }
    }

    flag=1;


	if (Ne <= 0) {
		fprintf(stderr, "Effective population size <= 0\n");
		flag = 0;
	}
	else fprintf(stderr, "Effective population size = %ld checked\n", Ne);

    if(THETA<0){
        fprintf(stderr, "Theta value < 0\n");
        flag=0;
    }
	else fprintf(stderr, "Theta value = %3.9lf checked\n", THETA);

    if(RHO<0){
        fprintf(stderr, "Rho value < 0\n");
        flag=0;
    }
	else fprintf(stderr, "Rho value = %3.9lf checked\n", RHO);

	if (SEQLEN < 0) {
		fprintf(stderr, "Sequence length < 0\n");
		flag = 0;
	}
	else fprintf(stderr, "Sequence length = %lf checked\n", SEQLEN);

	if (NumLineage < 0) {
		fprintf(stderr, "Number of lineages < 0\n");
		flag = 0;
	}
	else fprintf(stderr, "Number of lineages = %lf to stop\n", NumLineage);
    if(flag==0) abort();
}

void input2(FILE *infile)
{
    char string[200];
    char iden[2];
    short lines, flag, count;
    long i,j;
    long temp;
    double tot;
    
    L=0; /*number of loci*/
    M=0; /*number of genes*/
	Mu = 0;
    flag=1;
    lines=0;
    count=0;
	ntype_copy = 0;
    
    fprintf(stderr,"\n Data in file 2: \n");
    
    while(flag){
		if (fgets(string, sizeof(string), infile) == NULL) {
			fprintf(stderr, "Error in file \n");
			abort();
		}
        else{
            i=0;
            while (string[i]==' ' || string[i]=='\t') i++;
            iden[0]=string[i];
            switch(iden[0]){
                case EOF:
                    flag=0;
                    break;
                case '#':
                    flag=0;
                    break;
                case '\n':
                    break;

                case 'L':
                    if(L>0){
                        fprintf(stderr,"number of loci redefined on line %d. Ignored \n",lines);
                        break;
                    }
                    else{
                        while(string[i]!='=') i++;
                        i++;
                        sscanf_s(string+i,"%ld",&L);
                        fprintf(stderr,"number of loci %ld \n",L);
                        break;
                    }

                case 'G':
                    if(M>0){
                        fprintf(stderr,"number of genes redefined on line %d. Ignored \n",lines);
                        break;
                    }
                    else{
                        while(string[i]!='=') i++;
                        i++;
                        sscanf_s(string+i,"%ld",&M);
                        fprintf(stderr,"number of genes %ld \n",M);
                        break;
                    }

                case 'D':
                    if(ntype_copy>0){
                        fprintf(stderr,"number of distinct haplotypes on line %d. Ignored \n",lines);
                        break;
                    }
                    else{
                        while(string[i]!='=') i++;
                        i++;
                        sscanf_s(string+i,"%ld",&ntype_copy);
                        fprintf(stderr,"number of distinct haplotypes %ld \n",ntype_copy);
                        break;
                    }

				case 'S':
					if (Mu > 0) {
						fprintf(stderr, "Stationary distn redefined on line %d - ignored \n", lines);
						break;
					}
					else {
						while (string[i] != '=') i++;
						i++;
						sscanf_s(string + i, "%lf", &Mu);
						fprintf(stderr, "Stationary distn of allele 0 is %lf \n", Mu);
						break;
					}

                case 'H':
                    if(L==0 || ntype_copy==0){
                        fprintf(stderr,"Haplotypes redefined on line %d. Ignored \n",lines);
                        break;
                    }
                    else{
                        if((count&0x1)==0x1){
                            fprintf(stderr,"Haplotypes redefined on line %d. Ignored \n",lines);
                            break;
                        }
                        TYPE_copy=long_vec_init((L+1)*ntype_copy);
                        fprintf(stderr, "Haplotypes:\n");
                        for(j=0;j<(L+1)*ntype_copy;j++){
                            if(fscanf_s(infile,"%ld",&TYPE_copy[j])!=1){
                                fprintf(stderr,"Error in haplotypes for %ld th entry \n",j+1);
                                abort();
                            }
                            else{
                                if(j%(L+1)!=L) fprintf(stderr, "%ld ", TYPE_copy[j]);
                                else fprintf(stderr, "%ld\n", TYPE_copy[j]);
                            }
                        }
                        fprintf(stderr, "\n");
                        count+=1;
                        break;
                    }

                case 'M':
                    if((count&0x2)==0x2){
                        fprintf(stderr,"Mutation matrix redefined on line %d. Ignore \n",lines);
                        break;
                    }
                    fprintf(stderr, "Mutation matrix: ");
                    for(i=0;i<K;i++){
                        for(j=0;j<K;j++){
                            if(fscanf_s(infile,"%lf",&P[i][j])!=1){
                                fprintf(stderr,"Error in inputing mutation matrix element (%ld,%ld) \n",i+1,j+1);
                                abort();
                            }
                            else fprintf(stderr,"%lf ",P[i][j]);
                        }
                    }
                    fprintf(stderr, "\n");
                    count+=2;
                    break;

                case 'N':
                if(L==0){
                        fprintf(stderr,"New sequence redefined on line %d. Ignored \n",lines);
                        break;
                    }
                else{
                    if((count&0x4)==0x4){
                        fprintf(stderr,"New sequence redefined on line %d - ignored \n",lines);
                        break;
                    }
                    new = long_vec_init(L);
                    fprintf(stderr, "New sequence:");
                    for(i=0;i<L;i++){
                        if(fscanf_s(infile,"%ld",&new[i])!=1){
                            fprintf(stderr,"Error in new sequence for %ld th entry \n",i+1);
                            abort();
                        }
                        else fprintf(stderr, "%ld ", new[i]);
                    }
                    fprintf(stderr, "\n");
                    count+=4;
                    break;

                } 
                

                case 'P':
                    if(L==0){
                        fprintf(stderr, "Genetic positions defined before number of loci on line %d - ignore \n", lines);
                        break;
                    }
                    if((count&0x8)==0x8){
                        fprintf(stderr, "Genetic positions redefined on line %d \n", lines);
                        break;
                    }
					positions = dou_vec_init(L);
                    for(i=0;i<L;i++){
                        if(fscanf_s(infile, "%lf", &positions[i])!=1){
                            fprintf(stderr, "Error in inputting genetic positions on locus %ld \n", i+1);
                            abort();
                        }
                    }
                    count+=8;
                    break;

                default:
                    fprintf(stderr,"Non-standard input line %d in file 2. Ignored \n",lines);
                    break;
            }
            lines++;
        }
    }
    flag=1;
    if(L<=1){
        fprintf(stderr,"Number of loci <= 1 or undefined \n");
        flag=0;
    }
    if(ntype_copy <= 0){
        fprintf(stderr,"Number of haplotypes <= 0 or undefined \n");
        flag=0;
    }
    if(M <= 0){
        fprintf(stderr,"Number of genes <= 0 or undefined \n");
        flag=0;
    }
	if (Mu <= 0) {
		fprintf(stderr, "Stationary distn <= 0 or undefined \n");
		flag = 0;
	}
    if(flag==0) abort();

    /*check empirical distn */
    /*if((count&0x4)!=0x4){
        fprintf(stderr, "Empirical distn for N_tau not defined \n");
        flag=0;
    }
    else{
        for(i=0;i<2;i++){
            if(empirical[i]<=0){
                fprintf(stderr, "Mean and std of empirical distn should be positive \n");
                flag=0;
            }
        }
    }*/
    /*check positions */
    if((count&0x8)!=0x8){
        fprintf(stderr, "Genetic positions not defined \n");
        flag=0;
    }
    else{
        for(i=0;i<L;i++){
            if(positions[i]<0){
                fprintf(stderr, "Genetic position %ld should be non-negative \n", i+1);
                flag=0;
            }
        }
    }

    if(flag==0) abort();
    
    /*check TYPE*/
    temp=0;
    for(i=0;i<ntype_copy;i++){
        for(j=0;j<(L);j++){
            if(TYPE_copy[(L+1)*i+j]<0 || TYPE_copy[(L+1)*i+j]>K-1){
                fprintf(stderr,"Type for haplotype %ld at site %ld is not an allele \n",j+1,i+1);
                flag=0;
            }
        }
        if(TYPE_copy[(L+1)*i+L]<=0){
            fprintf(stderr,"Number of %ld th haplotype is less than or equal to 0 \n", i);
            flag=0;
        }
        temp+=TYPE_copy[(L+1)*i+L];
    }
    if(temp!=M){
        fprintf(stderr,"Number of genes from haplotypesdoes not equal the number of genes specified \n");
        flag=0;
    }

    /*check mutation matrix*/	
    if((count&0x2)!=0x2){
        fprintf(stderr,"Mutation matrix unspecified \n");
        flag=0;
    }
    else{
        for(i=0;i<K;i++){
            tot=0.0;
            for(j=0;j<K;j++){
                if(P[i][j]<0 || P[i][j]>1){
                    fprintf(stderr,"(%ld,%ld)th entry of mutation matrix not a probability\n",i+1,j+1);
                    flag=0;
                }
                tot+=P[i][j];
            }
            if(tot>=1.0000001 || tot<=0.999999) fprintf(stderr,"%ldth row of mutation matrix does not sum to one (renormalising)\n",i+1);
            for(j=0;j<K;j++) P[i][j]/=tot;
        }
    }
	

    if(flag==1) fprintf(stderr, "All checks done\n");
    if(flag==0) abort();
}
