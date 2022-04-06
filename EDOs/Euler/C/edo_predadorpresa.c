#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

int main(int argc, char* argv[]){
    int t=0, i=0, j=0;
    double deltaT = 0.;
    double T = 0; //Tfinal 
    double Ti = 0; 
    double passosDeTempo = 10000;
    
    /*FILE* entryfile;
    if ((entryfile = fopen(argv[1],"r")) == NULL){
        printf("Error on file opening! \n");
        exit(1);
    }*/
    
    printf("Digite os valores do T inicial, do T final e o número de passos de tempo:\n");
    scanf("%lf %lf %lf", &Ti, &T, &passosDeTempo);
    deltaT = (double)(T - Ti)/passosDeTempo;
    printf("deltaT = %lf\n", deltaT); 
    printf("Digite o número de iterações gravadas:\n");
    int numIteracoesGravadas = 100;    
    scanf("%d", &numIteracoesGravadas);
    int interval = (int)passosDeTempo/numIteracoesGravadas;  
    printf("intervalo : %d\n", interval);

    //Alocação e inicialização das populações 
    //Vetor de tamanho 2 para guardar os valores no passo de tempo anterior e no passo atual 
    double *h = (double*)malloc(2*sizeof(double)); //presa 
    double *p = (double*)malloc(2*sizeof(double)); //predador 
    for (int i = 0; i < 2; i++){
        h[i] = 0.;
        p[i] = 0.;
    }

    //Escreve o arquivo com os valores no tempo para os quais os resultados foram gravados 
    double time = 0.0;
    FILE* tempo = fopen("t.dat", "w");
    for (t=0; t <= passosDeTempo; t++){
        if ((t % interval) == 0)        
            fprintf(tempo, "%.2lf \n", time);        
        time = time + (double)deltaT;
    }
    fclose(tempo);

    FILE* htempo = fopen("presa.dat", "w");
    FILE* ptempo = fopen("predador.dat", "w");

    int t_ant = 0;
    int t_prox = 1;        

    //Condição inicial     
    h[t_ant] = 33;
    p[t_ant] = 6;
    //Parâmetros 
    double r = 0.09, kh = 1000, a = 0, m = 0.627; //r=0.937

    fprintf(htempo, "%lf\n", h[t_ant]); //Grava h[0]
    fprintf(ptempo, "%lf\n", p[t_ant]); //Grava p[0]

    //Loop temporal 
    for(int step = 1; step <= passosDeTempo; step++){
        //Troca os valores de t_prox e t_ant 
        if(step % 2 == 0) {
            t_prox = 0;
            t_ant = 1;
        }
        else {
            t_prox = 1;
            t_ant = 0;
        }
        //EDO presa 
        h[t_prox] = h[t_ant] + (r*h[t_ant] - a*h[t_ant]*p[t_ant])*deltaT;

        //EDO Predador
        p[t_prox] = p[t_ant] + ( a*h[t_ant]*p[t_ant] - m*p[t_ant])*deltaT;
        
        if((step % interval) == 0){
            fprintf(htempo, "%lf\n", h[t_prox]);
            fprintf(ptempo, "%lf\n", p[t_prox]);
        }
    }
    
    free(h);
    free(p);
    
    fclose(htempo);
    fclose(ptempo);
    //fclose(entryfile);
    return 0;
}
