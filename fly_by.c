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
// ....................................................................................................................
//      Constantes da simulação:
//  → Definições gerais.
#define N_DIMS 2                                        //  O número de dimensões do problema não deve ser alterado!!
                                                        //  Como o problema do Fly-By foi desenvolvido em cima do plano
                                                        //  de coordenadas polares, mudar esse valor não vai resultar
                                                        //  em alterações na simulação; porém, caso isso seja menor do
                                                        //  que 2, é provável que ocorra um erro na execução.

#define NUMERO_DE_TESTES 60                             //  Número de testes balísticos que serão realizados na simulação

#define STEPS_PARA_OUTPUT 900                           //  Número de "passos de integração" para que o programa exporte
                                                        //  um registro de saída. Sendo, por exemplo, dt = 0,2s, isso
                                                        //  resultaria num ponto a cada 180s (3 minutos) no arquivo
                                                        //  de saída.

#define DELTA_TIME 0.2                                  //  Step de integração numérica. (Em segundos)

#define TIME_FINAL 6.912e5                              //  Tempo máximo de integração. (Em segundos)
                                                        //  O valor de 6.912e5 segundos dá aproximadamente 8 dias.

//  → Constantes e massas.
#define CONSTANTE_GRAVITACIONAL 6.6743e-11              //  Constante gravitacional de Newton no sistema internacional.
#define MASSA_SOL 1.9885e30                             //  Massa do Sol. (Em quilogramas)
#define MASSA_MARTE 6.4171e23                           //  Massa de Marte. (Em quilogramas)

//  → Condições iniciais fixas
#define DISTANCIA_MARTE_SOL 2.2794e11                   //  Distância radial entre Marte e o Sol. (Em metros)
#define OFFSET_DA_SONDA 5e7                             //  Deslocamente entre a distância radial entre Marte e o Sol e
                                                        //  a posição mínima e máxima radial da sonda. (Em metros)

//  → Definições matemáticas
#define DEG_TO_RAD 0.0174532925                         //  Relação para converter graus para radianos.
#define RAD_TO_DEG 57.2957795                           //  Relação para converter radianos para graus.
// ....................................................................................................................
//      Funções auxiliares (aqui temos um "mini-header" dentro do arquivo *.c)
//  - Faz a simulação física do problema para as condições passadas nos argumentos.
//  const int test                          → Número de identificação do teste.
//  const double pos_ship_init              → Posição inicial no eixo 'x' da sonda.
//  double* d_min_value                     → Referência: altera a distância mínima encontrada entre a sonda e Marte.
//  double* b_min_value                     → Referência: altera o parâmetro de impacto mínimo encontrado.
//  double* delta_v_value                   → Referência: altera a variação de velocidade entre a sonda e Marte na entrada e saída.
//  double* deflection_angle_value          → Referência: altera o ângulo de deflexão encontrado.
//
//  * Os valores passados como "referência" são saídas da função simulate alterados em espaços de memória pré-alocados.
void simulate(int test, double pos_ship_init, double* d_min_value, double* b_min_value, double* delta_v_value, double* deflection_angle_value);
// ....................................................................................................................
//      Alocação global de memória:
char test_name[100];                                //  Nome da pasta onde os dados temporais serão salvos.

double theta_mars_init = 0.0;                       //  [º]          - Posição angular inicial de Marte. Esse valor
                                                    //              deve ser pequeno, algo entre (-0,25; 0,25).
double velocity_ship_init = 0.0;                    // [m/s]        - Velocidade radial inicial da sonda.
// ....................................................................................................................
//      Função de entrada do programa:
//  Compilação: gcc fly_by.c -lm -o fly_by
//  Permissão: chmod +x fly_by
//  Execução: ./fly_by [test_name] [mars_init_theta] [ship_init_radius_velocity]
int main(const int argc, const char *argv[]) {
    // ................................................................................................................
    //      Verifica o número de argumentos do input.
    if (argc != 4) {
        printf("Use: %s <test_name> <mars_init_theta> <ship_init_radius_velocity>\n",
            argv[0]);
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
    double pos_ship_init[NUMERO_DE_TESTES];             // [m]          - Posição inicial da sonda (sobre o eixo x)

    int i;                                              //              - Variável para iterações em primeiro nível.
    int j;                                              //              - Variável para iterações em segundo nível.

    double min_pos_ship_init;                           // [m]          - Valor mínimo para a posição radial inicial da sonda.
    double max_pos_ship_init;                           // [m]          - Valor máximo para a posição radial inicial da sonda.
    double pos_init_step;                               // [m]          - "Passo" entre os valores max e min.
    double progress;                                    //              - Contador de progresso.

    char filename[200];                                 //              - Nome do arquivo onde os dados globais serão salvos.
    FILE *fo;                                           //              - Ponteiro para o arquivo onde os dados serão salvos.
    // ................................................................................................................
    //      Salva o nome do teste numa variável global. Isso vai ser usado para o nome da pasta dos dados temporais,
    //  e também para o nome do arquivo de dados globais =D
    sprintf(test_name, "%s", argv[1]);

    //      Pega a posição angular de marte. E dá um aviso caso o valor seja muito grande.
    theta_mars_init = strtod(argv[2], NULL);
    if (fabs(theta_mars_init) > 0.25)
        printf("O valor do ângulo inicial de Marte está muito grande, isso pode resultar em uma interação fraca entre os corpos.\nTente usar algo entre (-0,25; 0,25).\n");
    theta_mars_init *= DEG_TO_RAD;

    //      Pega a velocidade radial inicial da sonda.
    velocity_ship_init = strtod(argv[3], NULL);

    //      Imprime um registro básico sobre o que o programa está a fazer.
    printf("Rodando o teste...\n");
    printf("Condições iniciais definidas pelo usuário: \n");
    printf("\t Distância radial entre Marte e o Sol: %.4e m (fixo por #define)\n", DISTANCIA_MARTE_SOL);
    printf("\t Posição angular inicial de Marte: %.6lf rad\n", theta_mars_init);
    printf("\t Velocidade radial inicial da sonda: %.2lf m/s\n", velocity_ship_init);
    // ................................................................................................................
    //      Calcula as posições iniciais dos testes.
    printf("Calculando posições iniciais para os testes...\n");
    printf("\t Serão realizados %d testes.\n", NUMERO_DE_TESTES);
    printf("\t Deslocamento de valores definido: %.4e m (fixo por #define)\n", OFFSET_DA_SONDA);

    min_pos_ship_init = DISTANCIA_MARTE_SOL - OFFSET_DA_SONDA;
    max_pos_ship_init = DISTANCIA_MARTE_SOL + OFFSET_DA_SONDA;
    printf("\t Serão adotados valores no intervalo: [%.4e m; %.4e m]\n", min_pos_ship_init, max_pos_ship_init);

    pos_init_step = (max_pos_ship_init - min_pos_ship_init) / (NUMERO_DE_TESTES - 1);
    printf("\t O passo de adição entre os valores assumido será de: %.8e m\n", pos_init_step);

    for (i = 0; i < NUMERO_DE_TESTES; i++) pos_ship_init[i] = min_pos_ship_init + (pos_init_step * i);
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
    //      Chama a função responsável pelas simulações numéricas de cada teste
    printf("\nRealizando simulações ... \n");
    for (i = 0; i < NUMERO_DE_TESTES; i++) {
        progress = (1.0 * i / NUMERO_DE_TESTES) * 100;

        simulate(i, pos_ship_init[i],
        &d_values[i], &b_values[i], &var_velocidade[i], &deflection_angle[i]);

        //      Barrinha de progresso =D
        printf("\r[");
        for (j = 0; j < NUMERO_DE_TESTES; j++) printf(j < i ? "#" : " ");
        printf("] %.2lf%%", progress);

        fflush(stdout);     //  Força a impressão =V
    }

    printf("\r[");
    for (i = 0; i < NUMERO_DE_TESTES; i++) printf("#");
    printf("] 100.00%%");
    fflush(stdout);
    printf("\n\n");
    // ................................................................................................................
    //      Salva os dados globais.
    printf("Salvando os dados globais em: '%s/global.csv'\n", test_name);

    sprintf(filename, "%s/global.csv", test_name);
    fo = fopen(filename, "w");

    //  - Cabeçalho do arquivo CSV.
    fprintf(fo, "i,pos_init,d_min,b_min,delta_v,deflection_angle\n");
    for (i = 0; i < NUMERO_DE_TESTES; i++) {
        fprintf(fo, "%d,%.15e,%.15e,%.15e,%.15e,%.15e\n",
            i + 1, pos_ship_init[i], d_values[i], b_values[i], var_velocidade[i], deflection_angle[i]);
    }

    fclose(fo);
    // ................................................................................................................
    printf("Simulação concluída =D\n\n");
    // ................................................................................................................
    return 0;
}
// ....................................................................................................................
//      Faz a simulação física para as condições iniciais apresentadas.
//  * Eu separei isso numa função a parte porque fica mais estruturado.
//  Aqui é o seguinte:
//  const int test                          → Número de identificação do teste.
//  const double pos_ship_init              → Posição inicial no eixo 'x' da sonda.
//  double* d_min_value                     → Referência: altera a distância mínima encontrada entre a sonda e Marte.
//  double* b_min_value                     → Referência: altera o parâmetro de impacto mínimo encontrado.
//  double* delta_v_value                   → Referência: altera a variação de velocidade entre a sonda e Marte na entrada e saída.
//  double* deflection_angle_value          → Referência: altera o ângulo de deflexão encontrado.
//
//  * Os valores passados como "referência" são saídas da função simulate alterados em espaços de memória pré-alocados.
void simulate(
    const int test,
    const double pos_ship_init,
    double* d_min_value,
    double* b_min_value,
    double* delta_v_value,
    double* deflection_angle_value) {
    // ................................................................................................................
    //          Declaração das variáveis locais.
    int f;                                              //                  - Contador para as saídas.

    //      Explicação do que farei abaixo:
    //  A integração numérica é em coordenadas polares, e é essencial eu utilizar elas aqui; porém, boa parte
    //  da análise numérica é feita com coordenadas cartesianas. Então é necesário calcular elas também.
    //  Por isso estarei a utilizar os dos sistemas de coordenadas.
    //
    double mars_coord_polar[N_DIMS + 1];                // [m, rad]         - Posição em coordenadas polares de Marte.
    double mars_coord_cartesian[N_DIMS + 1];            // [m, m]           - Posição em coordenadas cartesianas de Marte.
    double ship_coord_polar[N_DIMS + 1];                // [m, rad]         - Posição em coordenadas polares da sonda.
    double ship_coord_cartesian[N_DIMS + 1];            // [m, m]           - Posição em coordenadas cartesianas da sonda.
    double mars_velocity_polar[N_DIMS + 1];             // [m/s, rad/s]     - Velocidade em coordenadas polares de Marte.
    double mars_velocity_cartesian[N_DIMS + 1];         // [m/s, rad/s]     - Velocidade em coordenadas cartesianas de Marte.
    double ship_velocity_polar[N_DIMS + 1];             // [m/s, rad/s]     - Velocidade em coordenadas polares da sonda.
    double ship_velocity_cartesian[N_DIMS + 1];         // [m/s, rad/s]     - Velocidade em coordenadas cartesianas da sonda.

    //  Infelizmente eu também preciso de variáveis temporárias para armazenar as atualizações de estado...
    double mars_coord_polar_updated[N_DIMS + 1];        // [m, rad]         - Posição temporária atualizada pelo integrador de Marte.
    double ship_coord_polar_updated[N_DIMS + 1];        // [m, rad]         - Posição temporária atualizada pelo integrador da sonda.
    double ship_velocity_polar_updated[N_DIMS + 1];     // [m/s, rad/s]     - Velocidade temporária atualizada pelo integrador da sonda.

    //      Aqui temos algumas variáveis auxiliares para o cálculo das velocidades relativas (in, out)
    double relative_velocity_in[N_DIMS + 1];            // [m/s, m/s]       - Velocidade relativa de entrada (em coordenadas cartesianas)
    double relative_velocity_out[N_DIMS + 1];           // [m/s, m/s]       - Velocidade relativa de saída (em coordenadas cartesianas)

    double div;                                         // idk              - Variável auxiliar para o cálculo do divisor comum nas Eqs. de integração.
    double time;                                        // [s]              - Tempo de intregrassão.
    double relative_distance;                           // [m]              - Distância relativa entre a sonda e Marte.
    double impact_parameter;                            // [m]              - Parâmetro de impacto entre a sonda e Marte.
    char filename[200];                                 //                  - Arquivos onde os dados da simulação serão salvos.
    FILE *fo;                                           //                  - Ponteiro de acesso para o arquivo.
    // ................................................................................................................
    //          Condições iniciais, aplicadas...
    mars_coord_polar[1] = DISTANCIA_MARTE_SOL;          //  Posição radial inicial de Marte (isso é fixo!)
    mars_coord_polar[2] = theta_mars_init;              //  Posição angular inicial de Marte
    mars_velocity_polar[1] = 0.0;                       //  Velocidade racial incial de Marte (sempre zero!)
    mars_velocity_polar[2]                              //  Velocidade angular inicial de Marte (fixo também!)
        = sqrt(CONSTANTE_GRAVITACIONAL * MASSA_SOL / (DISTANCIA_MARTE_SOL * DISTANCIA_MARTE_SOL * DISTANCIA_MARTE_SOL));
    ship_coord_polar[1] = pos_ship_init;                //  Posição radial inicial da sonda (definido na entrada da função)
    ship_coord_polar[2] = 0.0;                          //  Posição angular inicial da sonda (fixo também!)
    ship_velocity_polar[1] = velocity_ship_init;        //  Velocidade radial inicial da sonda
    ship_velocity_polar[2]                              //  Velocidade angular inicial da sonda (a fórmula é padrão, mas no caso depende da posição radial inicial)
        = sqrt(CONSTANTE_GRAVITACIONAL * MASSA_SOL / (pos_ship_init * pos_ship_init * pos_ship_init));

    //      Converte as coordenadas cartesianas iniciais.
    //  - Posição de Marte (x e y, em ordem)
    mars_coord_cartesian[1] = mars_coord_polar[1] * cos(mars_coord_polar[2]);
    mars_coord_cartesian[2] = mars_coord_polar[1] * sin(mars_coord_polar[2]);
    //  - Posição da sonda (x e y, em ordem)
    ship_coord_cartesian[1] = ship_coord_polar[1] * cos(ship_coord_polar[2]);
    ship_coord_cartesian[2] = ship_coord_polar[1] * sin(ship_coord_polar[2]);
    //  - Velocidade de Marte (x e y, em ordem)
    mars_velocity_cartesian[1] = mars_velocity_polar[1] * cos(mars_coord_polar[2]) - mars_velocity_polar[2] * sin(mars_coord_polar[2]);
    mars_velocity_cartesian[2] = mars_velocity_polar[1] * sin(mars_coord_polar[2]) + mars_velocity_polar[2] * cos(mars_coord_polar[2]);
    //  - Velocidade da sonda (x e y, em ordem)
    ship_velocity_cartesian[1] = ship_velocity_polar[1] * cos(ship_coord_polar[2]) - ship_velocity_polar[2] * sin(ship_coord_polar[2]);
    ship_velocity_cartesian[2] = ship_velocity_polar[1] * sin(ship_coord_polar[2]) + ship_velocity_polar[2] * cos(ship_coord_polar[2]);

    //      Calcula a distância relativa e o parâmetro de impacto iniciais.
    //  Isso é uma conta gigante, mas representa o valor de |\vec d| apresentado no trabalho escrito; e também da Eq.(33).
    relative_distance = sqrt(ship_coord_polar[1] * ship_coord_polar[1] + mars_coord_polar[1] * mars_coord_polar[1] - 2 * ship_coord_polar[1] * mars_coord_polar[1] * cos(ship_coord_polar[2] - mars_coord_polar[2]));
    impact_parameter = ((ship_coord_cartesian[1] - mars_coord_cartesian[1]) * (ship_velocity_cartesian[2] - mars_velocity_cartesian[2]) - (ship_coord_cartesian[2] - mars_coord_cartesian[2]) * (ship_velocity_cartesian[1] - mars_velocity_cartesian[1])) /
        sqrt((ship_velocity_cartesian[1] - mars_velocity_cartesian[1]) * (ship_velocity_cartesian[1] - mars_velocity_cartesian[1]) + (ship_velocity_cartesian[2] - mars_velocity_cartesian[2]) * (ship_velocity_cartesian[2] - mars_velocity_cartesian[2]));

    *d_min_value = relative_distance;
    *b_min_value = impact_parameter;

    //      Calcula a velocidade relativa de entrada
    relative_velocity_in[1] = ship_velocity_cartesian[1] - mars_velocity_cartesian[1];
    relative_velocity_in[2] = ship_velocity_cartesian[2] - mars_velocity_cartesian[2];
    // ................................................................................................................
    //          Prepara para salvar os dados.
    f = 0;
    sprintf(filename, "%s/data_%03d.csv", test_name, test + 1);
    fo = fopen(filename, "w");
    //  Como eu não vou usar GnuPlot, vou adicionar um cabeçalho no arquivo; e meio que o formato acaba virando um CSV.
    fprintf(fo, "t,x_mars,y_mars,x_ship,y_ship,v_x_mars,v_y_mars,v_x_ship,v_y_ship,d,b\n");
    // ................................................................................................................
    //          Processo de simulação numérica.
    for (time = 0.0; time < TIME_FINAL; time += DELTA_TIME) { // NOLINT(*-flp30-c)
        // ............................................................................................................
        //          Adiciona os dados ao arquivo de saída.
        if (f <= 0) {
            f = STEPS_PARA_OUTPUT;
            fprintf(fo, "%.8e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e,%.15e\n",
                time, mars_coord_cartesian[1], mars_coord_cartesian[2], ship_coord_cartesian[1], ship_coord_cartesian[2],
                mars_velocity_cartesian[1], mars_velocity_cartesian[2], ship_velocity_cartesian[1], ship_velocity_cartesian[2],
                relative_distance, impact_parameter);
        }

        f--;
        // ............................................................................................................
        //          Realiza a integração numérica em coordenadas polares, conforme Eqs~(26-30).
        //  - Calcula o divisor comum das expressões. Essa expressão é a mesma usada no cálculo de 'relative_distance'.
        //  A expressão é: div = r_s^2 + r_M^2 - 2 * r_s * r_M * cos(\theta_s - \theta_M)
        div = ship_coord_polar[1] * ship_coord_polar[1] + mars_coord_polar[1] * mars_coord_polar[1] - 2 * ship_coord_polar[1] * mars_coord_polar[1] * cos(ship_coord_polar[2] - mars_coord_polar[2]);
        //  - Ajusta a posição angular de Marte.
        mars_coord_polar_updated[2] = mars_coord_polar[2] + mars_velocity_polar[2] * DELTA_TIME;
        //  - Ajusta a posição da sonda.
        ship_coord_polar_updated[1] = ship_coord_polar[1] + ship_velocity_polar[1] * DELTA_TIME;
        ship_coord_polar_updated[2] = ship_coord_polar[2] + ship_velocity_polar[2] * DELTA_TIME;
        //  - Ajusta a velocidade da sonda.
        ship_velocity_polar_updated[1] = ship_velocity_polar[1] + (ship_coord_polar[1] * ship_velocity_polar[2] * ship_velocity_polar[2] -
            CONSTANTE_GRAVITACIONAL * MASSA_SOL / (ship_coord_polar[1] * ship_coord_polar[1]) +
            CONSTANTE_GRAVITACIONAL * MASSA_MARTE * (ship_coord_polar[1] - mars_coord_polar[1] * cos(ship_coord_polar[2] - mars_coord_polar[2])) /
            (div * sqrt(div))) * DELTA_TIME;
        ship_velocity_polar_updated[2] = ship_velocity_polar[2] +
            (CONSTANTE_GRAVITACIONAL * MASSA_MARTE * mars_coord_polar[1] * sin(ship_coord_polar[2] - mars_coord_polar[2]) /
                (ship_coord_polar[1] * div * sqrt(div)) - 2 * ship_velocity_polar[1] * ship_velocity_polar[2] / ship_coord_polar[1]) * DELTA_TIME;
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
        mars_velocity_cartesian[1] = mars_velocity_polar[1] * cos(mars_coord_polar[2]) - mars_velocity_polar[2] * sin(mars_coord_polar[2]);
        mars_velocity_cartesian[2] = mars_velocity_polar[1] * sin(mars_coord_polar[2]) + mars_velocity_polar[2] * cos(mars_coord_polar[2]);
        //  - Velocidade da sonda (x e y, em ordem)
        ship_velocity_cartesian[1] = ship_velocity_polar[1] * cos(ship_coord_polar[2]) - ship_velocity_polar[2] * sin(ship_coord_polar[2]);
        ship_velocity_cartesian[2] = ship_velocity_polar[1] * sin(ship_coord_polar[2]) + ship_velocity_polar[2] * cos(ship_coord_polar[2]);

        //      Calcula a distância relativa e o parâmetro de impacto iniciais.
        //  Isso é uma conta gigante, mas representa o valor de |\vec d| apresentado no trabalho escrito; e também da Eq.(33).
        relative_distance = sqrt(ship_coord_polar[1] * ship_coord_polar[1] + mars_coord_polar[1] * mars_coord_polar[1] - 2 * ship_coord_polar[1] * mars_coord_polar[1] * cos(ship_coord_polar[2] - mars_coord_polar[2]));
        impact_parameter = ((ship_coord_cartesian[1] - mars_coord_cartesian[1]) * (ship_velocity_cartesian[2] - mars_velocity_cartesian[2]) - (ship_coord_cartesian[2] - mars_coord_cartesian[2]) * (ship_velocity_cartesian[1] - mars_velocity_cartesian[1])) /
            sqrt((ship_velocity_cartesian[1] - mars_velocity_cartesian[1]) * (ship_velocity_cartesian[1] - mars_velocity_cartesian[1]) + (ship_velocity_cartesian[2] - mars_velocity_cartesian[2]) * (ship_velocity_cartesian[2] - mars_velocity_cartesian[2]));

        if (relative_distance < *d_min_value) *d_min_value = relative_distance;
        if (impact_parameter < *b_min_value) *b_min_value = impact_parameter;
        // ............................................................................................................
    }
    // ................................................................................................................
    //          Fecha o arquivo de dados.
    fclose(fo);
    // ................................................................................................................
    //          Calcula o ângulo de deflexão e a variação da velocidade relativa.
    //      Primeiro, vou calcular a velocidade relativa de saída...
    relative_velocity_out[1] = ship_velocity_cartesian[1] - mars_velocity_cartesian[1];
    relative_velocity_out[2] = ship_velocity_cartesian[2] - mars_velocity_cartesian[2];

    //      Agora calculamos a variação da velocidade relativa, ou seja, a diferença dos módulos de in e out.
    *delta_v_value = sqrt(relative_velocity_out[1] * relative_velocity_out[1] + relative_velocity_out[2] * relative_velocity_out[2]) -
        sqrt(relative_velocity_in[1] * relative_velocity_in[1] + relative_velocity_in[2] * relative_velocity_in[2]);

    //      E agora, calculamos o ângulo de deflexão...
    *deflection_angle_value = (relative_velocity_in[1] * relative_velocity_out[1] + relative_velocity_in[2] * relative_velocity_out[2]) /
        (sqrt(relative_velocity_in[1] * relative_velocity_in[1] + relative_velocity_in[2] * relative_velocity_in[2]) *
            sqrt(relative_velocity_out[1] * relative_velocity_out[1] + relative_velocity_out[2] * relative_velocity_out[2]));
    if (*deflection_angle_value > 1.0) *deflection_angle_value = 1.0;
    if (*deflection_angle_value < -1.0) *deflection_angle_value = -1.0;
    *deflection_angle_value = acos(*deflection_angle_value);
    // ................................................................................................................
}