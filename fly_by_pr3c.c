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
                                                        //  Como o problema do Fly-By foi desenvolvido em cima do plano
                                                        //  de coordenadas polares, mudar esse valor não vai resultar
                                                        //  em alterações na simulação; porém, caso isso seja menor do
                                                        //  que 2, é provável que ocorra um erro na execução.

#define NUMERO_DE_TESTES 240                            //  Número de testes balísticos que serão realizados na simulação

#define STEPS_PARA_OUTPUT 180                           //  Valor base para definir a cada quantos steps de integração
                                                        //  ocorrerá a exportação de dados e checagem dos critérios de parada.
                                                        //  Pense que para dt = 1 segundo isso é feito a cada STEPS_PARA_OUTPUT segundos.

//  → Constantes e massas.
#define CONSTANTE_GRAVITACIONAL 6.6743e-11              //  Constante gravitacional de Newton no sistema internacional.
#define MASSA_SOL 1.9885e30                             //  Massa do Sol. (Em quilogramas)
#define MASSA_MARTE 6.4171e23                           //  Massa de Marte. (Em quilogramas)

//  → Condições iniciais fixas
#define DISTANCIA_MARTE_SOL 2.2794e11                   //  Distância radial entre Marte e o Sol. (Em metros)
#define RAIO_MARTE 3.3895E6                             //  Raio do planeta Marte. (Em metros)

//  → Definições matemáticas
#define DEG_TO_RAD 0.0174532925                         //  Relação para converter graus para radianos.
#define RAD_TO_DEG 57.2957795                           //  Relação para converter radianos para graus.
// ....................................................................................................................
//      Funções auxiliares (aqui temos um "mini-header" dentro do arquivo *.c)
//  - Faz a simulação física do problema para as condições passadas nos argumentos.
//  const int test                          → Número de identificação do teste.
//  const double b                          → Parâmetro de impacto.
//  double* d_min_value                     → Referência: altera a distância mínima encontrada entre a sonda e Marte.
//  double* delta_v_value                   → Referência: altera a variação de velocidade heliocêntrica da sonda.
//  double* delta_v_value_rel               → Referência: altera a variação de velocidade relativa da sonda.
//  double* deflection_angle_value          → Referência: altera o ângulo de deflexão encontrado.
//  int* collision                          → Referência: altera o indicador de colisão. É 0 caso não tenha ocorrido colisão; e 1 caso contrário.
//  double* time_end                        → Referência: altera o tempo em segundo que levou para finalizar a simulação.
//
//  * Os valores passados como "referência" são saídas da função simulate alterados em espaços de memória pré-alocados.
void simulate(int test, double b, double* d_min_value, double* delta_v_value, double* delta_v_value_rel, double* deflection_angle_value, int* collision, double* time_end);

//  - Isso aqui é só uma função extra para converter o tempo de ETA de segundos para um formato "melhor"
//  Pode ignorar =D
void format_time(double seconds, char *buffer);
// ....................................................................................................................
//      Alocação global de memória:
char test_name[100];                                //  Nome da pasta onde os dados temporais serão salvos.

double mars_angle_init;                             //  Posição angular inicial de Marte.
double r_factor;                                    //  Fator que multiplica R_Marte para definir a esfera de influência do planeta.
double v_sonda_init;                                //  Módulo da velocidade inicial da sonda.
double max_int_time;                                //  Tempo total de simulação (critério de parada de emergência)
double dt;                                          //  Timestep de integração.
double stop_value;                                  //  Fator de parada da simulação. (em relação a órbitas de Marte)

int steps_to_output;                                //  Passos de integração para a exportação.
// ....................................................................................................................
//      Função de entrada do programa:
//  Compilação: gcc fly_by.c -lm -o fly_by
//  Permissão: chmod +x fly_by
//  Execução: ./fly_by [test_name] [x_init_factor] [mars_init_angle] [velocity_infinity] [b_min_factor] [b_max_factor] [max_time] [dt]
int main(const int argc, const char *argv[]) {
    // ................................................................................................................
    //      Verifica o número de argumentos do input.
    if (argc != 9) {
        printf("Use: %s <test_name> <x_init_factor> <velocity_infinity> <b_min_factor> <b_max_factor> <max_time> <dt>\n",
            argv[0]);
        printf("- <test_name>: Nome do teste, e da pasta onde a saída será salva.\n");
        printf("- <r_factor>: Fator multiplicando R_Marte que define o raio de influência do planeta. No trabalho usamos um valor igual a 50.\n");
        printf("- <mars_init_angle>: Ângulo inicial de Marte no sistema de coordenadas cartesiano referenciado no Sol, em graus. No trabalho usamos um valor igual a -0.01\n");
        printf("- <velocity_infinity>: Velocidade da sonda no infinito, em metros por segundo. No trabalho usamos o valor de 2600 m/s.\n");
        printf("- <b_min_factor>: Fator mínimo usado na definição do intervalo de valores para o parâmetro de impacto. No trabalho usamos um valor igual a -10\n");
        printf("- <b_max_factor>: Fator máximo usado na definição do intervalo de valores para o parâmetro de impacto. No trabalho usamos um valor igual a 10.\n");
        printf("- <max_time>: Critério de parada de emergência. É o tempo máximo que pode ser gasto com a integração antes dela ser abortada, sem segundos. No trabalho foi utilizado 10e10 segundos.\n");
        printf("- <dt>: Passo temporal utilizado na integração, em segundos. Não deve ser muito grande já que é usado o método de Euler. No trabalho foi utilizado 0,001 s.\n");
        return 1;
    }
    // ................................................................................................................
    //      Declara as variáveis locais.
    double d_values[NUMERO_DE_TESTES];                  // [m]          - Distância relativa mínima entre a sonda e
                                                        //              Marte registrada para um teste.
    double b_values[NUMERO_DE_TESTES];                  // [m]          - Parâmetro de impacto usado no teste.
    double var_velocidade_rel[NUMERO_DE_TESTES];        // [m/s]        - Variação do módulo da velocidade relativa de
                                                        //              entrada e de saída da sonda no processo.
    double var_velocidade_helio[NUMERO_DE_TESTES];      // [m/s]        - Variação do módulo da velocidade heliocêntrica da sonda.
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
    r_factor = strtod(argv[2], NULL);
    r_factor *= RAIO_MARTE;
    stop_value = r_factor;

    //      Pega a inclinação inicial de Marte (em radianos)
    mars_angle_init = strtod(argv[3], NULL);
    mars_angle_init *= DEG_TO_RAD;

    //      Pega a velocidade da sonda no infinito, e calcula a velocidade inicial da sonda.
    v_sonda_init = strtod(argv[4], NULL);
    v_sonda_init = sqrt(v_sonda_init * v_sonda_init + 2 * CONSTANTE_GRAVITACIONAL * MASSA_MARTE / fabs(r_factor));

    //      Calcula o intervalo de valores de parâmetro de impacto que serão utilizados.
    min_b_factor = strtod(argv[5], NULL);
    max_b_factor = strtod(argv[6], NULL);

    min_b_factor *= RAIO_MARTE;
    max_b_factor *= RAIO_MARTE;

    if (fabs(r_factor) < max_b_factor) {
        printf("O raio de influência da esfera não pode ser menor do que o fator de impacto máximo.\n");
        return 1;
    }

    b_step = (max_b_factor - min_b_factor) / (NUMERO_DE_TESTES - 1);

    //      Tempo e passo de integração.
    max_int_time = strtod(argv[7], NULL);
    dt = strtod(argv[8], NULL);

    steps_to_output = (int) (STEPS_PARA_OUTPUT / dt);

    //      Calcula os valores de fator de impacto que serão usados.
    for (i = 0; i < NUMERO_DE_TESTES; i++) b_values[i] = min_b_factor + b_step * i;

    // ................................................................................................................
    //      Imprime um registro básico sobre o que o programa está a fazer.
    printf("Rodando o teste...\n");
    printf("Condições iniciais definidas: \n");
    printf("\t Massa do Sol utilizada: %.4e kg \n", MASSA_SOL);
    printf("\t Raio de Marte utilizado: %.4e metros \n", RAIO_MARTE);
    printf("\t Massa de Marte utilizada: %.4e kg \n", MASSA_MARTE);
    printf("\t Raio da órbita de Marte utilizada: %.4e metros\n", DISTANCIA_MARTE_SOL);
    printf("\t Valor raio de influência da esfera é de: %.4e metros\n", r_factor);
    printf("\t Valor do parâmetro de impacto pertencente ao intervalo [%.4e m; %.4e m], com passo igual a %.4e metros\n", min_b_factor, max_b_factor, b_step);
    printf("\t Valor do módulo da velocidade inicial da sonda: %.4e metros por segundo\n", v_sonda_init);
    printf("\t Tempo máximo de integração: %.4e segundos\n", max_int_time);
    printf("\t Passo de integração: %.4lf s\n", dt);
    // ................................................................................................................
    //      Cria a pasta onde as coisas serão salvas.
    //  !! Esse trecho do código funciona apenas no MacOS e no Linux. Isso não é aplicável no Windows.
    //  No Windows é necessário substituir essa implementação com o uso da biblioteca 'direct.h'.
    sprintf(filename, "%s/pr3c", test_name);
    if (mkdir(filename, 0755) == 0) printf("Pasta do problema de 3 corpos: '%s'\n", filename);
    else {
        perror("Falha ao criar o diretório para os arquivos do problema de 3 corpos. Verifique se a pasta já existe, caso isso seja verdade, delete-a ou renomei-a.");
        return 1;
    }
    // ................................................................................................................
    //      Chama a função responsável pelas simulações numéricas de cada teste.
    printf("\nRealizando simulações ... \n");
    begin = clock();
    for (i = 0; i < NUMERO_DE_TESTES; i++) {
        progress = (1.0 * i / NUMERO_DE_TESTES) * 100;

        //      Chama a simulação para o parâmetro de impacto b_values[i].
        simulate(i, b_values[i], &d_values[i], &var_velocidade_helio[i], &var_velocidade_rel[i], &deflection_angle[i], &collision[i], &times[i]);

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

    sprintf(filename, "%s/global_pr3c.csv", test_name);
    fo = fopen(filename, "w");

    //  - Cabeçalho do arquivo CSV.
    fprintf(fo, "i,b,d_min,delta_v,delta_v_rel,deflection_angle,collision,t\n");
    for (i = 0; i < NUMERO_DE_TESTES; i++) {
        fprintf(fo, "%d,%.15e,%.15e,%.15e,%.15e,%.15e,%d,%.15e\n",
            i + 1, b_values[i], d_values[i], var_velocidade_helio[i], var_velocidade_rel[i], deflection_angle[i] * RAD_TO_DEG, collision[i],times[i]);
    }

    fclose(fo);
    // ................................................................................................................
    printf("Simulação concluída =D\n\n");
    // ................................................................................................................
    return 0;
}
// ....................................................................................................................
//  - Faz a simulação física do problema para as condições passadas nos argumentos.
//  const int test                          → Número de identificação do teste.
//  const double b                          → Parâmetro de impacto.
//  double* d_min_value                     → Referência: altera a distância mínima encontrada entre a sonda e Marte.
//  double* delta_v_value                   → Referência: altera a variação de velocidade heliocêntrica da sonda.
//  double* delta_v_value_rel               → Referência: altera a variação de velocidade relativa da sonda.
//  double* deflection_angle_value          → Referência: altera o ângulo de deflexão encontrado.
//  int* collision                          → Referência: altera o indicador de colisão. É 0 caso não tenha ocorrido colisão; e 1 caso contrário.
//  double* time_end                        → Referência: altera o tempo em segundo que levou para finalizar a simulação.
//
//  * Os valores passados como "referência" são saídas da função simulate alterados em espaços de memória pré-alocados.
void simulate(int test, double b, double* d_min_value, double* delta_v_value, double* delta_v_value_rel, double* deflection_angle_value, int* collision, double* time_end) {
    // ................................................................................................................
    //          Declaração das variáveis locais.
    int f;                                              //                  - Contador para as saídas.

    //      Variáveis de estado (coordenadas polares)
    double mars_coord_polar[N_DIMS + 1];                // [m, rad]         - Posição em coordenadas polares de Marte.
    double mars_velocity_polar[N_DIMS + 1];             // [m/s, rad/s]     - Velocidade em coordenadas polares de Marte.
    double ship_coord_polar[N_DIMS + 1];                // [m, rad]         - Posição em coordenadas polares da sonda.
    double ship_velocity_polar[N_DIMS + 1];             // [m/s, rad/s]     - Velocidade em coordenadas polares da sonda.

    //      Variáveis de estado (coordenadas cartesianas)
    double mars_coord_cartesian[N_DIMS + 1];            // [m, m]           - Posição em coordenadas cartesianas de Marte.
    double mars_velocity_cartesian[N_DIMS + 1];         // [m/s, m/s]       - Velocidade em coordenadas cartesianas de Marte.
    double ship_coord_cartesian[N_DIMS + 1];            // [m, m]           - Posição em coordenadas cartesianas da sonda.
    double ship_velocity_cartesian[N_DIMS + 1];         // [m/s, m/s]       - Velocidade em coordenadas cartesianas da sonda.

    //  Também preciso de variáveis temporárias para armazenar as atualizações de estado...
    double mars_coord_polar_updated[N_DIMS + 1];        // [m, rad]         - Posição temporária atualizada pelo integrador de Marte.
    double ship_coord_polar_updated[N_DIMS + 1];        // [m, rad]         - Posição temporária atualizada pelo integrador da sonda.
    double ship_velocity_polar_updated[N_DIMS + 1];     // [m/s, rad/s]     - Velocidade temporária atualizada pelo integrador da sonda.

    //      Vetores de velocidade de entrada e saída.
    double velocity_in[N_DIMS + 1];                     // [m/s, m/s]       - Vetor de velocidade de entrada. (referencial do Sol)
    double velocity_out[N_DIMS + 1];                    // [m/s, m/s]       - Vetor de velocidade de saída. (referecial do Sol)
    double velocity_in_rel[N_DIMS + 1];                 // [m/s, m/s]       - Vetor de velocidade de entrada. (referencial de Marte)
    double velocity_out_rel[N_DIMS + 1];                // [m/s, m/s]       - Vetor de velocidade de saída. (referecial de Marte)
    double velocity_out_direction[N_DIMS + 1];          // [m/s, m/s]       - Vetor unitário de direção do vetor de saída no ponto de parada. Isso define a direção do velocity_out.

    double div;                                         // [mˆ2]            - Fator comum de divisão. É o módulo quadrado do vetor distância.
    double time;                                        // [s]              - Tempo de intregrassão.
    double distance;                                    // [m]              - Distância relativa entre a sonda e Marte.
    char filename[200];                                 //                  - Arquivos onde os dados da simulação serão salvos.
    FILE *fo;                                           //                  - Ponteiro de acesso para o arquivo.
    // ................................................................................................................
    //          Condições iniciais, aplicadas...
    //  1. Começando por marte, o raio da órbita do planeta é tabelado; e o ângulo foi fornecido.
    mars_coord_polar[1] = DISTANCIA_MARTE_SOL;
    mars_coord_polar[2] = mars_angle_init;
    mars_velocity_polar[1] = 0.0;
    mars_velocity_polar[2] = sqrt(CONSTANTE_GRAVITACIONAL * MASSA_SOL / (DISTANCIA_MARTE_SOL * DISTANCIA_MARTE_SOL * DISTANCIA_MARTE_SOL));

    //  2. Converte as coordenadas de marte para os valores cartesianos.
    mars_coord_cartesian[1] = mars_coord_polar[1] * cos(mars_coord_polar[2]);
    mars_coord_cartesian[2] = mars_coord_polar[1] * sin(mars_coord_polar[2]);
    mars_velocity_cartesian[1] = - mars_coord_polar[1] * mars_velocity_polar[2] * sin(mars_coord_polar[2]);
    mars_velocity_cartesian[2] = mars_coord_polar[1] * mars_velocity_polar[2] * cos(mars_coord_polar[2]);

    //  3. Calcula a posição cartesiana inicial da sonda.
    ship_coord_cartesian[1] = mars_coord_cartesian[1] + sqrt(r_factor * r_factor - b * b) * sin(mars_coord_polar[2]) - b * cos(mars_coord_polar[2]);
    ship_coord_cartesian[2] = mars_coord_cartesian[2] - sqrt(r_factor * r_factor - b * b) * cos(mars_coord_polar[2]) - b * sin(mars_coord_polar[2]);

    //  4. Converte os valores de posição para coordendas polares.
    ship_coord_polar[1] = sqrt(ship_coord_cartesian[1] * ship_coord_cartesian[1] + ship_coord_cartesian[2] * ship_coord_cartesian[2]);
    ship_coord_polar[2] = atan2(ship_coord_cartesian[2], ship_coord_cartesian[1]);

    //  5. Calcula a velocidade relativa de entrada; e aproveita para calcular o valor da velocidade heliocêntrica.
    velocity_in_rel[1] = - v_sonda_init * sin(mars_coord_polar[2]);
    velocity_in_rel[2] = v_sonda_init * cos(mars_coord_polar[2]);

    ship_velocity_cartesian[1] = velocity_in_rel[1] + mars_velocity_cartesian[1];
    ship_velocity_cartesian[2] = velocity_in_rel[2] + mars_velocity_cartesian[2];

    velocity_in[1] = ship_velocity_cartesian[1];
    velocity_in[2] = ship_velocity_cartesian[2];

    //  6. Converte as velocidades cartesianas para polar.
    ship_velocity_polar[1] = (ship_coord_cartesian[1] * ship_velocity_cartesian[1] + ship_coord_cartesian[2] * ship_velocity_cartesian[2]) / ship_coord_polar[1];
    ship_velocity_polar[2] = (ship_coord_cartesian[1] * ship_velocity_cartesian[2] - ship_coord_cartesian[1] * ship_velocity_cartesian[1]) / (ship_coord_polar[1] * ship_coord_polar[1]);

    //  7. Distância entre a sonda e Marte.
    distance = sqrt((ship_coord_cartesian[1] - mars_coord_cartesian[1]) * (ship_coord_cartesian[1] - mars_coord_cartesian[1]) + (ship_coord_cartesian[2] - mars_coord_cartesian[2]) * (ship_coord_cartesian[2] - mars_coord_cartesian[2]));

    //  8. Configurações adicionais...
    *collision = 0;
    *d_min_value = distance;
    // ................................................................................................................
    //          Prepara para salvar os dados.
    f = 0;
    sprintf(filename, "%s/pr3c/data_%03d.csv", test_name, test + 1);
    fo = fopen(filename, "w");
    //  Como eu não vou usar GnuPlot, vou adicionar um cabeçalho no arquivo; e meio que o formato acaba virando um CSV.
    fprintf(fo, "t,x_mars,y_mars,x_ship,y_ship,v_x_mars,v_y_mars,v_x_ship,v_y_ship,d\n");
    // ................................................................................................................
    //          Processo de simulação numérica.
    for (time = 0; time < max_int_time; time += dt) { // NOLINT(*-flp30-c)
        // ............................................................................................................
        //          Adiciona os dados ao arquivo de saída.
        if (f <= 0) {
            f = steps_to_output;
            fprintf(fo, "%.8e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                time, mars_coord_cartesian[1], mars_coord_cartesian[2], ship_coord_cartesian[1], ship_coord_cartesian[2], mars_velocity_cartesian[1], mars_velocity_cartesian[2], ship_velocity_cartesian[1], ship_velocity_cartesian[2], distance);

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
        //          Realiza a integração numérica, conforme Eqs~(34-38).
        div = ship_coord_polar[1] * ship_coord_polar[1] + mars_coord_polar[1] * mars_coord_polar[1] - 2 * ship_coord_polar[1] * mars_coord_polar[1] * cos(ship_coord_polar[2] - mars_coord_polar[2]);

        mars_coord_polar_updated[2] = mars_coord_polar[2] + mars_velocity_polar[2] * dt;
        ship_coord_polar_updated[1] = ship_coord_polar[1] + ship_velocity_polar[1] * dt;
        ship_coord_polar_updated[2] = ship_coord_polar[2] + ship_velocity_polar[2] * dt;
        ship_velocity_polar_updated[1] = ship_velocity_polar[1] + (ship_coord_polar[1] * ship_velocity_polar[2] * ship_velocity_polar[2] -
            CONSTANTE_GRAVITACIONAL * MASSA_SOL / (ship_coord_polar[1] * ship_coord_polar[1]) -
            CONSTANTE_GRAVITACIONAL * MASSA_MARTE * (ship_coord_polar[1] - mars_coord_polar[1] * cos(ship_coord_polar[2] - mars_coord_polar[2])) /
            (div * sqrt(div))) * dt;
        ship_velocity_polar_updated[2] = ship_velocity_polar[2] +
            (- CONSTANTE_GRAVITACIONAL * MASSA_MARTE * mars_coord_polar[1] * sin(ship_coord_polar[2] - mars_coord_polar[2]) /
                (ship_coord_polar[1] * div * sqrt(div)) - 2 * ship_velocity_polar[1] * ship_velocity_polar[2] / ship_coord_polar[1]) * dt;

        //  - Aplica as atualizações de posição e velocidade.
        mars_coord_polar[2] = mars_coord_polar_updated[2];
        ship_coord_polar[1] = ship_coord_polar_updated[1];
        ship_coord_polar[2] = ship_coord_polar_updated[2];
        ship_velocity_polar[1] = ship_velocity_polar_updated[1];
        ship_velocity_polar[2] = ship_velocity_polar_updated[2];
        // ............................................................................................................
        //          Atualiza os dados de coordenadas cartesianas, etc.
        //      Converte as coordenadas cartesianas iniciais.
        //  - Posição de Marte (x e y, em ordem)
        mars_coord_cartesian[1] = mars_coord_polar[1] * cos(mars_coord_polar[2]);
        mars_coord_cartesian[2] = mars_coord_polar[1] * sin(mars_coord_polar[2]);
        //  - Posição da sonda (x e y, em ordem)
        ship_coord_cartesian[1] = ship_coord_polar[1] * cos(ship_coord_polar[2]);
        ship_coord_cartesian[2] = ship_coord_polar[1] * sin(ship_coord_polar[2]);
        //  - Velocidade de Marte (x e y, em ordem)
        mars_velocity_cartesian[1] = mars_velocity_polar[1] * cos(mars_coord_polar[2]) - mars_coord_polar[1] * mars_velocity_polar[2] * sin(mars_coord_polar[2]);
        mars_velocity_cartesian[2] = mars_velocity_polar[1] * sin(mars_coord_polar[2]) + mars_coord_polar[1] * mars_velocity_polar[2] * cos(mars_coord_polar[2]);
        //  - Velocidade da sonda (x e y, em ordem)
        ship_velocity_cartesian[1] = ship_velocity_polar[1] * cos(ship_coord_polar[2]) - ship_coord_polar[1] * ship_velocity_polar[2] * sin(ship_coord_polar[2]);
        ship_velocity_cartesian[2] = ship_velocity_polar[1] * sin(ship_coord_polar[2]) + ship_coord_polar[1] * ship_velocity_polar[2] * cos(ship_coord_polar[2]);
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
    //          Calcula o ângulo de deflexão e a variação da velocidade.
    //  →   Vamos começar calculando a variação do módulo da velocidade heliocêntrica. Para isso vamos pegar a
    //  velocidade heliocêntrica de saída:
    velocity_out[1] = ship_velocity_cartesian[1];
    velocity_out[2] = ship_velocity_cartesian[2];
    //  Agora comparamos os módulos...
    *delta_v_value = sqrt(velocity_out[1] * velocity_out[1] + velocity_out[2] * velocity_out[2]) - sqrt(velocity_in[1] * velocity_in[1] + velocity_in[2] * velocity_in[2]);

    //  →   Na sequência, calculamos a variação da velocidade relativa; é a mesma coisa, mas a gente precisa descontar a velocidade
    //  de marte...
    velocity_out_rel[1] = ship_velocity_cartesian[1] - mars_velocity_cartesian[1];
    velocity_out_rel[2] = ship_velocity_cartesian[2] - mars_velocity_cartesian[2];

    *delta_v_value_rel = sqrt(velocity_out_rel[1] * velocity_out_rel[1] + velocity_out_rel[2] * velocity_out_rel[2]) - sqrt(velocity_in_rel[1] * velocity_in_rel[1] + velocity_in_rel[2] * velocity_in_rel[2]);

    //  →   Vamos agora pegar a direção do vetor velocidade... (aqui vou usar cálculo repetido, mas é para deixar mais claro oq estou tentando fazer)
    div = sqrt(velocity_out_rel[1] * velocity_out_rel[1] + velocity_out_rel[2] * velocity_out_rel[2]);
    velocity_out_direction[1] = velocity_out[1] / div;
    velocity_out_direction[2] = velocity_out[2] / div;

    //  → E agora, calculamos o ângulo de deflexão...
    *deflection_angle_value = (velocity_in_rel[1] * velocity_out_rel[1] + velocity_in_rel[2] * velocity_out_rel[2]) /
        (sqrt(velocity_in_rel[1] * velocity_in_rel[1] + velocity_in_rel[2] * velocity_in_rel[2]) *
            sqrt(velocity_out_rel[1] * velocity_out_rel[1] + velocity_out_rel[2] * velocity_out_rel[2]));

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