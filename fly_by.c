//
//  File in Portuguese.
//  Author: Gabriel V. Ferreira (2025)
//  GitHub: https://github.com/gabriel-ferr/Fly-by
// ....................................................................................................................
//      Comentários para desabilitar algumas funções de análise do CLion:
// ReSharper disable CppJoinDeclarationAndAssignment
// ....................................................................................................................
//      Bibliotecas:
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
// ....................................................................................................................
//      Constantes da simulação:
//  → Definições gerais.
#define N_DIMS 2                                        //  O número de dimensões do problema não deve ser alterado!!
                                                        //  O jeito em que o código foi desenvolvido é em cima de um
                                                        //  plano bidimensional cartesiano; aumentar isso não vai resultar
                                                        //  em nada, além de alocação extra de memória. Porém, diminuir pode
                                                        //  resultar em erros de acesso de memória.

#define NUMERO_DE_TESTES 240                            //  Número de testes balísticos que serão realizados na simulação

#define STEPS_PARA_OUTPUT 180                           //  Valor base para definir a cada quantos steps de integração
                                                        //  ocorrerá a exportação de dados e checagem dos critérios de parada.
                                                        //  Pense que para dt = 1 segundo isso é feito a cada STEPS_PARA_OUTPUT segundos.

//  → Constantes e massas.
#define CONSTANTE_GRAVITACIONAL 6.6743e-11              //  Constante gravitacional de Newton no sistema internacional.
#define MASSA_MARTE 6.4171e23                           //  Massa de Marte. (Em quilogramas)

//  → Condições iniciais fixas
#define RAIO_MARTE 3.3895E6                             //  Raio do planeta Marte. (Em metros)

//  → Matemática
#define RAD_TO_DEG 57.2957795                           //  Converte radianos para graus.
// ....................................................................................................................
//      Funções auxiliares (aqui temos um "mini-header" dentro do arquivo *.c)
//  - Faz a simulação física do problema para as condições passadas nos argumentos.
//  const int test                          → Número de identificação do teste.
//  const double b                          → Parâmetro de impacto.
//  double* d_min_value                     → Referência: altera a distância mínima encontrada entre a sonda e Marte.
//  double* delta_v_value                   → Referência: altera a variação de velocidade entre a sonda e Marte na entrada e saída.
//  double* deflection_angle_value          → Referência: altera o ângulo de deflexão encontrado.
//  int* collision                          → Referência: altera o indicador de colisão. É 0 caso não tenha ocorrido colisão; e 1 caso contrário.
//  double* time_end                        → Referência: altera o tempo em segundo que levou para finalizar a simulação.
//
//  * Os valores passados como "referência" são saídas da função simulate alterados em espaços de memória pré-alocados.
void simulate(int test, double b, double* d_min_value, double* delta_v_value, double* deflection_angle_value, int* collision, double* time_end);

//  - Isso aqui é só uma função extra para converter o tempo de ETA de segundos para um formato "melhor"
//  Pode ignorar =D
void format_time(double seconds, char *buffer);
// ....................................................................................................................
//      Alocação global de memória:
char test_name[100];                                //  Nome da pasta onde os dados temporais serão salvos.

double x_init;                                      //  Valor inicial no eixo x.
double v_infinite_in;                               //  Velocidade da sonda no infinito.
double v_x_init;                                    //  Velocidade da sonda inicial (apenas no eixo x)
double max_int_time;                                //  Tempo total de simulação (critério de parada de emergência)
double dt;                                          //  Timestep de integração.
double stop_value;                                  //  Fator de parada da simulação. (em relação a órbitas de Marte)

int steps_to_output;                                //  Passos de integração para a exportação.
// ....................................................................................................................
//      Função de entrada do programa:
//  Compilação: gcc fly_by.c -lm -o fly_by
//  Permissão: chmod +x fly_by
//  Execução: ./fly_by [test_name] [x_init_factor] [velocity_infinity] [b_min_factor] [b_max_factor] [max_time] [dt]
int main(const int argc, const char *argv[]) {
    // ................................................................................................................
    //      Verifica o número de argumentos do input.
    if (argc != 8) {
        printf("Use: %s <test_name> <x_init_factor> <velocity_infinity> <b_min_factor> <b_max_factor> <max_time> <dt>\n",
            argv[0]);
        printf("- <test_name>: Nome do teste, e da pasta onde a saída será salva.\n");
        printf("- <x_init_factor>: Fator multiplicando R_Marte na definição de x(0). No trabalho usamos um valor igual a 20. Também é usado no critério de parada.\n");
        printf("- <velocity_infinity>: Velocidade da sonda no infinito, em metros por segundo. No trabalho usamos o valor de 3000 m/s.\n");
        printf("- <b_min_factor>: Fator mínimo usado na definição do intervalo de valores para o parâmetro de impacto. No trabalho usamos um valor igual a 2. Não recomendo tomar um valor menor do que 2.\n");
        printf("- <b_max_factor>: Fator máximo usado na definição do intervalo de valores para o parâmetro de impacto. No trabalho usamos um valor igual a 6. Não recomendo tomar um valor maior do que 10.\n");
        printf("- <max_time>: Critério de parada de emergência. É o tempo máximo que pode ser gasto com a integração antes dela ser abortada, sem segundos. No trabalho foi utilizado 1,8e5 segundos.\n");
        printf("- <dt>: Passo temporal utilizado na integração, em segundos. Não deve ser muito grande já que é usado o método de Euler. No trabalho foi utilizado 0,05 segundos.\n");
        return 1;
    }
    // ................................................................................................................
    //      Declara as variáveis locais.
    double d_values[NUMERO_DE_TESTES];                  // [m]          - Distância relativa mínima entre a sonda e
                                                        //              Marte registrada para um teste.
    double b_values[NUMERO_DE_TESTES];                  // [m]          - Parâmetro de impacto mínimo registrado para
                                                        //              um teste.
    double var_velocidade[NUMERO_DE_TESTES];            // [m/s]        - Módulo da variação da velocidade relativa de
                                                        //              entrada e de saída da sonda no processo.
    double deflection_angle[NUMERO_DE_TESTES];          // [º]          - Ângulo de deflexão entre a sonda e Marte.
    double times[NUMERO_DE_TESTES];                     // [s]          - Tempo total que cada simulação utilizou, em segundos.
    int collision[NUMERO_DE_TESTES];                    //              - Indicador de colisão. É 1 caso tenha ocorrido
                                                        //              - uma colisão da sonda com Marte, e 0 caso contrário.

    int i;                                              //              - Variável para iterações em primeiro nível.
    int j;                                              //              - Variável para iterações em segundo nível.

    double min_b_factor;                                // [m]          - Fator mínimo para o parâmetro b
    double max_b_factor;                                // [m]          - Fator máximo para o parâmetro b
    double b_step;                                      // [m]          - "Passo" entre os valores max e min de b.
    double progress;                                    //              - Contador de progresso.
    double elapsed;                                     //              - Tempo decorrido da simulação.
    double total_time;                                  //              - Tempo estimado total de simulação.
    double remaining;                                   //              - Tempo restante de simulação.
    clock_t begin;                                      //              - Momento em que o processo começou.
    char elapsed_str[50];                               //              - “String” com o tempo atual de processamento.
    char remaining_str[50];                             //              - "String" com o tempo restante de processamento.

    char filename[200];                                 //              - Nome do arquivo onde os dados globais serão salvos.
    FILE *fo;                                           //              - Ponteiro para o arquivo onde os dados serão salvos.
    // ................................................................................................................
    //      Salva o nome do teste numa variável global. Isso vai ser usado para o nome da pasta dos dados temporais,
    //  e também para o nome do arquivo de dados globais =D
    sprintf(test_name, "%s", argv[1]);

    //      Cálcula o valor de x(0) a partir do fator passado na chamada do código.
    x_init = strtod(argv[2], NULL);
    x_init *= -RAIO_MARTE;
    stop_value = (-1) * x_init;        //  É a mesma coisa que o x_init, mas positivo (e não necessariamente sobre x)

    //      Pega a velocidade da sonda no infinito, e calcula a velocidade inicial da sonda.
    v_infinite_in = strtod(argv[3], NULL);
    v_x_init = sqrt(v_infinite_in * v_infinite_in + 2 * CONSTANTE_GRAVITACIONAL * MASSA_MARTE / fabs(x_init));

    //      Calcula o intervalo de valores de parâmetro de impacto que serão utilizados.
    min_b_factor = strtod(argv[4], NULL);
    max_b_factor = strtod(argv[5], NULL);

    if (min_b_factor > max_b_factor) {
        printf("O fator mínimo (argv[4]: %d) precisa ser menor do que o fator máximo (argv[5]: %d).\n", (int) min_b_factor, (int) max_b_factor);
        return 1;
    }

    min_b_factor *= RAIO_MARTE;
    max_b_factor *= RAIO_MARTE;

    if (fabs(x_init) < max_b_factor) {
        printf("A posição inicial não pode ser menor do que o fator de impacto máximo.\n");
        return 1;
    }

    b_step = (max_b_factor - min_b_factor) / (NUMERO_DE_TESTES - 1);

    //      Tempo e passo de integração.
    max_int_time = strtod(argv[6], NULL);
    dt = strtod(argv[7], NULL);

    steps_to_output = (int) (STEPS_PARA_OUTPUT / dt);

    //      Calcula os valores de fator de impacto que serão usados.
    for (i = 0; i < NUMERO_DE_TESTES; i++) b_values[i] = min_b_factor + b_step * i;
    // ................................................................................................................
    //      Imprime um registro básico sobre o que o programa está a fazer.
    printf("Rodando o teste...\n");
    printf("Condições iniciais definidas: \n");
    printf("\t Raio de Marte utilizado: %.4e metros \n", RAIO_MARTE);
    printf("\t Massa de Marte utilizada: %.4e metros \n", MASSA_MARTE);
    printf("\t Valor de x(0): %.4e metros\n", x_init);
    printf("\t Valor de y(0) pertencente ao intervalo [%.4e m; %.4e m], com passo igual a %.4e metros\n", min_b_factor, max_b_factor, b_step);
    printf("\t Valor de vx(0): %.4e metros por segundo\n", v_x_init);
    printf("\t Valor de vy(0): %.4e metros por segundo\n", 0.0);
    printf("\t Tempo máximo de integração: %.4e segundos\n", max_int_time);
    printf("\t Passo de integração: %.4lf s\n", dt);
    // ................................................................................................................
    //      Cria a pasta onde as coisas serão salvas.
    //  !! Esse trecho do código funciona apenas no MacOS e no Linux. Isso não é aplicável no Windows.
    //  No Windows é necessário substituir essa implementação com o uso da biblioteca 'direct.h'.
    if (mkdir(test_name, 0755) == 0) printf("Os dados serão salvos na pasta: '%s'\n", test_name);
    else {
        perror("Falha ao criar o diretório do teste. Verifique se a pasta já existe, caso isso seja verdade, delete-a ou renomei-a.");
        return 1;
    }
    // ................................................................................................................
    //      Chama a função responsável pelas simulações numéricas de cada teste.
    printf("\nRealizando simulações ... \n");
    begin = clock();
    for (i = 0; i < NUMERO_DE_TESTES; i++) {
        progress = (1.0 * i / NUMERO_DE_TESTES) * 100;

        //      Chama a simulação para o parâmetro de impacto b_values[i].
        simulate(i, b_values[i], &d_values[i], &var_velocidade[i], &deflection_angle[i], &collision[i], &times[i]);

        //      Barrinha de progresso (modo avançado com ETA) =D
        elapsed = (double) (clock() - begin) / CLOCKS_PER_SEC;
        total_time = (elapsed / (i > 0 ? i : 1)) * NUMERO_DE_TESTES;
        remaining = total_time - elapsed;

        format_time(elapsed, elapsed_str);
        format_time(remaining, remaining_str);

        printf("\r");       //  Limpa
        for (j = 0; j < 220; j++) printf(" ");
        fflush(stdout);

        printf("\r[");
        for (j = 0; j < 100; j++) printf(progress >= j ? "#" : " ");
        printf("] %.2lf%%, Elapsed: %s, ETA: %s", progress, elapsed_str, remaining_str);

        fflush(stdout);     //  Força a impressão =V
    }

    printf("\r");       //  Limpa
    for (j = 0; j < 220; j++) printf(" ");
    fflush(stdout);

    printf("\r[");
    for (i = 0; i < 100; i++) printf("#");
    printf("] 100.00%%, Total time: %s", elapsed_str);
    fflush(stdout);
    printf("\n\n");
    // ................................................................................................................
    //      Salva os dados globais.
    printf("Salvando os dados globais em: '%s/global.csv'\n", test_name);

    sprintf(filename, "%s/global.csv", test_name);
    fo = fopen(filename, "w");

    //  - Cabeçalho do arquivo CSV.
    fprintf(fo, "i,b,d_min,delta_v,deflection_angle,collision,t\n");
    for (i = 0; i < NUMERO_DE_TESTES; i++) {
        fprintf(fo, "%d,%.15e,%.15e,%.15e,%.15e,%d,%.15e\n",
            i + 1, b_values[i], d_values[i], var_velocidade[i], deflection_angle[i] * RAD_TO_DEG, collision[i],times[i]);
    }

    fclose(fo);
    // ................................................................................................................
    printf("Simulação concluída =D\n\n");
    // ................................................................................................................
    return 0;
}
// ....................................................................................................................
//  * Eu separei isso numa função a parte porque fica mais organizado.
//  - Faz a simulação física do problema para as condições passadas nos argumentos.
//  const int test                          → Número de identificação do teste.
//  const double b                          → Parâmetro de impacto.
//  double* d_min_value                     → Referência: altera a distância mínima encontrada entre a sonda e Marte.
//  double* delta_v_value                   → Referência: altera a variação de velocidade entre a sonda e Marte na entrada e saída.
//  double* deflection_angle_value          → Referência: altera o ângulo de deflexão encontrado.
//
//  * Os valores passados como "referência" são saídas da função simulate alterados em espaços de memória pré-alocados.
//  double* time_end                        → Referência: altera o tempo em segundo que levou para finalizar a simulação.
void simulate(const int test, const double b, double* d_min_value, double* delta_v_value, double* deflection_angle_value, int* collision, double* time_end) {
    // ................................................................................................................
    //          Declaração das variáveis locais.
    int f;                                              //                  - Contador para as saídas.

    //      Variáveis de estado.
    double r[N_DIMS + 1];                               // [m, m]           - Posição da sonda.
    double v[N_DIMS + 1];                               // [m/s, m/s]       - Velocidade da sonda.
    //  → A aceleração foi incorporada diretamente nas equações de movimento.
    //      Variáveis de estado temporárias.
    double r_temp[N_DIMS + 1];                          // [m, m]           - Posição temporária da sonda.
    double v_temp[N_DIMS + 1];                          // [m/s, m/s]       - Velocidade temporária da sonda.

    //      Vetores de velocidade de entrada e saída.
    double velocity_in[N_DIMS + 1];                     // [m/s, m/s]       - Vetor de velocidade relativa de entrada. (No infinito)
    double velocity_out[N_DIMS + 1];                    // [m/s, m/s]       - Vetor de velocidade relativa de saída. (No infinito, cálculado com correção de energia)
    double velocity_out_direction[N_DIMS + 1];          // [m/s, m/s]       - Vetor unitário de direção do vetor de saída no ponto de parada. Isso define a direção do velocity_out.

    double div;                                         // [mˆ2]            - Fator comum de divisão. É o módulo quadrado do vetor distância.
    double time;                                        // [s]              - Tempo de intregrassão.
    double distance;                                    // [m]              - Distância relativa entre a sonda e Marte.
    char filename[200];                                 //                  - Arquivos onde os dados da simulação serão salvos.
    FILE *fo;                                           //                  - Ponteiro de acesso para o arquivo.
    // ................................................................................................................
    //          Condições iniciais, aplicadas...
    r[1] = x_init;
    r[2] = b;
    v[1] = v_x_init;
    v[2] = 0.0;

    velocity_in[1] = v_infinite_in;
    velocity_in[2] = 0.0;

    distance = sqrt(r[1] * r[1] + r[2] * r[2]);

    *collision = 0;
    *d_min_value = distance;
    // ................................................................................................................
    //          Prepara para salvar os dados.
    f = 0;
    sprintf(filename, "%s/data_%03d.csv", test_name, test + 1);
    fo = fopen(filename, "w");
    //  Como eu não vou usar GnuPlot, vou adicionar um cabeçalho no arquivo; e meio que o formato acaba virando um CSV.
    fprintf(fo, "t,x,y,v_x,v_y,d\n");
    // ................................................................................................................
    //          Processo de simulação numérica.
    for (time = 0; time < max_int_time; time += dt) { // NOLINT(*-flp30-c)
        // ............................................................................................................
        //          Adiciona os dados ao arquivo de saída.
        if (f <= 0) {
            f = steps_to_output;
            fprintf(fo, "%.8e,%.12e,%.12e,%.12e,%.12e,%.12e\n",
                time, r[1], r[2], v[1], v[2], distance);

            //      Critérios de parada.
            //  1. Verifica se a sonda colidiu com Marte.
            if (distance < RAIO_MARTE) {
                *d_min_value = distance;
                *collision = 1;
                break;
            }

            //  2. Verifica se a sonda está suficientemente longe de Marte.
            //  O critério temporal aqui é apenas para impedir que ele pare no começo da simulação.
            if (distance >= stop_value && time > 10 * STEPS_PARA_OUTPUT) {
                break;
            }
        }

        f--;
        // ............................................................................................................
        //          Realiza a integração numérica, conforme Eqs~(13-16).
        div = r[1] * r[1] + r[2] * r[2];
        r_temp[1] = r[1] + v[1] * dt;
        r_temp[2] = r[2] + v[2] * dt;
        v_temp[1] = v[1] - (CONSTANTE_GRAVITACIONAL * MASSA_MARTE * r[1] / (div * sqrt(div))) * dt;
        v_temp[2] = v[2] - (CONSTANTE_GRAVITACIONAL * MASSA_MARTE * r[2] / (div * sqrt(div))) * dt;

        //          Atualiza os estados.
        r[1] = r_temp[1];
        r[2] = r_temp[2];
        v[1] = v_temp[1];
        v[2] = v_temp[2];
        // ............................................................................................................
        //          Calcula a distância entre Marte e a sonda.
        distance = sqrt(div);
        if (distance < *d_min_value) *d_min_value = distance;
        // ............................................................................................................
    }
    // ................................................................................................................
    //          Fecha o arquivo de dados.
    fclose(fo);
    // ................................................................................................................
    //          Calcula o ângulo de deflexão e a variação da velocidade relativa.
    //  → Começamos com a velocidade de saída no infinito.
    //  A velocidade atualmente armazenada em 'v' ainda apresenta um erro devido ao campo potenical de Marte; então,
    //  a gente vai corrigir esse valor para o infinito, descontando a parcela associada a energia potencial gravitacional
    //  com base na posição 'r' em que 'v' foi obtida.
    //  O primeiro passo aqui é pegar a direção de 'v':
    div = sqrt(v[1] * v[1] + v[2] * v[2]);
    velocity_out_direction[1] = v[1] / div;
    velocity_out_direction[2] = v[2] / div;
    //  Na sequência, a gente calcula o módulo corrigido da velocidade no infinito.
    div = sqrt(div * div - 2 * CONSTANTE_GRAVITACIONAL * MASSA_MARTE / sqrt(r[1] * r[1] + r[2] * r[2]));    // 'div' dentro da expressão é o módulo da velocidade.
    velocity_out[1] = div * velocity_out_direction[1];
    velocity_out[2] = div * velocity_out_direction[2];

    //  → Como temos os vetores de entrada e saída no infinito, calculamos a variação dos módulos...
    *delta_v_value = sqrt(velocity_out[1] * velocity_out[1] + velocity_out[2] * velocity_out[2]) - sqrt(velocity_in[1] * velocity_in[1] + velocity_in[2] * velocity_in[2]);

    //  → E agora, calculamos o ângulo de deflexão...
    *deflection_angle_value = (velocity_in[1] * velocity_out[1] + velocity_in[2] * velocity_out[2]) /
        (sqrt(velocity_in[1] * velocity_in[1] + velocity_in[2] * velocity_in[2]) *
            sqrt(velocity_out[1] * velocity_out[1] + velocity_out[2] * velocity_out[2]));

    if (*deflection_angle_value > 1) *deflection_angle_value = 1.0;
    if (*deflection_angle_value < -1) *deflection_angle_value = -1.0;

    *deflection_angle_value = acos(*deflection_angle_value);

    //  → Por fim, seta o tempo total usado para a integração.
    *time_end = time;
}
// ....................................................................................................................
//      ** Função para mostrar o tempo no ETA em segundos, minutos, etc.
//  Pode ignorar isso aqui; não dei muita atenção em manter organizado também...
void format_time(double seconds, char *buffer) {
    const int sec = (int) seconds;
    if (sec < 60) {
        sprintf(buffer, "%d segundos", sec);
    } else if (sec < 3600) {
        const int m = sec / 60;
        const int s = sec % 60;
        sprintf(buffer, "%d minutos e %d segundos", m, s);
    } else if (sec < 86400) {
        const int h = sec / 3600;
        const int m = sec % 3600;
        sprintf(buffer, "%d horas e %d minutos", h, m);
    } else {
        const int d = sec / 86400;
        const int h = sec % 86400;
        sprintf(buffer, "%d dias e %d minutos", d, h);
    }
}
// ....................................................................................................................