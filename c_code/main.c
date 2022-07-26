#include <stdio.h>
#include <stdlib.h>


float quimiotaxia(float ponto_atual, float valor_maximo){
    return ponto_atual/(valor_maximo + ponto_atual);
}

float gradiente(float ponto_posterior, float ponto_anterior, float valor_maximo){
    return (quimiotaxia(ponto_posterior, valor_maximo) - quimiotaxia(ponto_anterior, valor_maximo));
}

float f_func(float populacao, float valor_maximo){
    return populacao*populacao/(valor_maximo + populacao);
}


//Fazer struct com os campos desejados e a funcao init retorna uma instancia dessa struct. A funcao model vai operar em cima dessa struct
void init_model(){

}

void model(){

}

int main(){

    return 0;
}