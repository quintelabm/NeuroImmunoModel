#include <stdio.h>
#include <stdlib.h>
#define T_final (int)28
#define X_final (int)20
#define ht 0.0002
#define hx 0.2
int num_points_x = (int) X_final/hx;
int num_points_t = (int) T_final/ht;

double quimiotaxia(double ponto_atual, double valor_maximo){
    return ponto_atual/(valor_maximo + ponto_atual);
}

double gradiente(double ponto_posterior, double ponto_anterior, double valor_maximo){
    return (quimiotaxia(ponto_posterior, valor_maximo) - quimiotaxia(ponto_anterior, valor_maximo));
}

double f_func(double populacao, double valor_maximo){
    return populacao*populacao/(valor_maximo + populacao);
}

struct model_def {
    double** odc_new;
    double** odc_old;
    double** mic_new;
    double** mic_old;
    double** tke_new;
    double** tke_old;
    double** ant_new;
    double** ant_old;
    double** ddc_new;
    double** ddc_old;
    double** dda_new;
    double** dda_old;
    double num_figuras;
    double V_BV; //Volume(valor total) de blood vessels
    double V_LV; //Volume(valor total) de perivascular area
    

};

struct model_def mdl;

//Fazer struct com os campos desejados e a funcao init retorna uma instancia dessa struct. A funcao model vai operar em cima dessa struct
void initialize_model(){
    double** matrix_aux[num_points_x][num_points_x];
    mdl.odc_new = matrix_aux;
    mdl.odc_old = matrix_aux;
    mdl.mic_new = matrix_aux;
    mdl.mic_old = matrix_aux;
    mdl.tke_new = matrix_aux;
    mdl.tke_old = matrix_aux;
    mdl.ant_new = matrix_aux;
    mdl.ant_old = matrix_aux;
    mdl.ddc_new = matrix_aux;
    mdl.ddc_old = matrix_aux;
    mdl.dda_new = matrix_aux;
    mdl.dda_old = matrix_aux;
}

void model(){

}

int main(){

    return 0;
}