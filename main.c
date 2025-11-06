//
//  File in Portuguese.
//  Author: Gabriel V. Ferreira (2025)
//  --------------------------------------------------
#include <math.h>
#include <stdio.h>

//      Configurações constantes.
#define G_CONST 6.6743e-11                      //              -   Constante gravitacional de Newton.
#define SUN_M 1.9885e30                         //  [kg]        -   Massa do Sol.
#define MARS_M 6.4171e23                        //  [kg]        -   Massa de Marte.
#define MARS_R 2.28e11                          //  [m]         -   Raio da órbita de Marte com relação ao Sol na origem.
#define MARS_INIT_THETA 45                      //  [º]         -   Ângulo inicial de Marte, em graus.
#define MARS_INIT_V_THETA 1.058e-7              //  [rad/s]     -
#define SHIP_R 2.2799e11                        //  [m]         -   Raio da órbita da Sonda com relação ao Sol. Esse valor é dado pela
                                                //                  condicional 3 do trabalho.
#define SHIP_INIT_V_R 100                       //  [m/s]       -   Velocidade radial inicial da Sonda.
#define SHIP_INIT_V_THETA 1.058e-7              //  [rad/s]     -   Velocidade angular inicial da Sonda. Dado pela condicional 6 do trabalho.
#define HILL_R 1.0843e9                         //  [m]         -   Raio da esfera de Hill, para o critério de finalização da simulação.
#define MAX_INTEGRATION_TIME 1e6                //  [s]         -   Tempo máximo de integração da simulação.
#define DELTA_T 0.5                             //  [s]         -   Step de integração.
#define B_MIN_VALUE 8e6                         //  [m]         -   Valor mínimo para o parâmetro de impacto.
#define B_MAX_VALUE 2.1e9                       //  [m]         -   Valor máximo para o parâmetro de impacto.
#define NUMBER_OF_TESTS 100                     //              -   Número de valores de parâmetro de impacto que serão usados / testes realizados.
#define FILE_SAVE 600                           //              -   Número de steps de integração para salvar uma saída.

double module_v(const double v_s, const double theta_s, const double r_s, const double omega_s) {
    const double v_x = v_s * cos(theta_s) - r_s * omega_s * sin(theta_s);
    const double v_y = v_s * sin(theta_s) + r_s * omega_s * cos(theta_s);

    return sqrt(v_x * v_x + v_y * v_y);
}

double comp_deflection(const double theta_init_value, const double r_s_out, const double v_r_out, const double theta_s_out,
    const double omega_s_out, const double theta_M_out) {

    double mars_theta_in = MARS_INIT_THETA * (M_PI / 180.0);

    //      Velocidade da sonda com relação ao sol em coordenadas cartesianas.
    const double v_ship_x_in = SHIP_INIT_V_R * cos(theta_init_value) - SHIP_R * SHIP_INIT_V_THETA * sin(theta_init_value);
    const double v_ship_y_in = SHIP_INIT_V_R * sin(theta_init_value) + SHIP_R * SHIP_INIT_V_THETA * cos(theta_init_value);

    const double v_ship_x_out = v_r_out * cos(theta_s_out) - r_s_out * omega_s_out * sin(theta_s_out);
    const double v_ship_y_out = v_r_out * sin(theta_s_out) + r_s_out * omega_s_out * cos(theta_s_out);

    //      Velocidade de Marte com relação ao sol em coordenadas cartesianas.
    const double v_mars_x_in = - MARS_R * MARS_INIT_V_THETA * sin(mars_theta_in);
    const double v_mars_y_in = MARS_R * MARS_INIT_V_THETA * cos(mars_theta_in);

    const double v_mars_x_out = - MARS_R * MARS_INIT_V_THETA * sin(theta_M_out);
    const double v_mars_y_out = MARS_R * MARS_INIT_V_THETA * cos(theta_M_out);

    //      Velocidades relativas da sonda em relação a Marte.
    const double v_rel_in_x = v_ship_x_in - v_mars_x_in;
    const double v_rel_in_y = v_ship_y_in - v_mars_y_in;
    const double v_rel_out_x = v_ship_x_out - v_mars_x_out;
    const double v_rel_out_y = v_ship_y_out - v_mars_y_out;

    //      Produtos internos e normas.
    const double dot = v_rel_in_x * v_rel_out_x + v_rel_in_y * v_rel_out_y;
    const double mag_in = sqrt(v_rel_in_x * v_rel_in_x + v_rel_in_y * v_rel_in_y);
    const double mag_out = sqrt(v_rel_out_x * v_rel_out_x + v_rel_out_y * v_rel_out_y);

    //      Cálculo do ângulo.
    double cos_delta = dot / (mag_in * mag_out);
    if (cos_delta > 1.0) cos_delta = 1.0;
    if (cos_delta < -1.0) cos_delta = -1.0;

    double delta = acos(cos_delta);
    const double cross = v_rel_in_x * v_rel_out_y - v_rel_in_y * v_rel_out_x;

    if (cross < 0) delta = -delta;
    return delta;
}

void simulate(const int i, double *rel_velocity, double *deflection, double theta_init_value) {
    int f = 0;

    double mars_theta = MARS_INIT_THETA * (M_PI / 180.0);
    double ship_r = SHIP_R;
    double ship_v = SHIP_INIT_V_R;
    double ship_theta = theta_init_value;
    double ship_v_theta = SHIP_INIT_V_THETA;

    const double v_in = module_v(ship_v, ship_theta, ship_r, ship_v_theta);

    char filename[100];
    sprintf(filename, "../data/result_%d.dat", i);          //      Note que isso está adaptado para o CMake!!
    FILE *fo = fopen(filename, "w");

    for (double t = 0.0; t < MAX_INTEGRATION_TIME; t += DELTA_T) {
        if (f <= 0) {
            fprintf(fo, "%e %e %e %e %e %e\n", t, mars_theta, ship_r, ship_theta, ship_v, ship_v_theta);
            f = FILE_SAVE;

            //      Testa o critério de saída.
            const double r_rel_mag = sqrt(ship_r * ship_r + MARS_R * MARS_R - 2 * ship_r * MARS_R * cos(ship_theta - mars_theta));
            if (r_rel_mag > 2.0 * HILL_R) break;
        }

        f--;

        //      Dividor comum.
        const double div = (ship_r * ship_r) + (MARS_R * MARS_R) - 2 * ship_r * MARS_R * cos(ship_theta - mars_theta);

        //      Posição angular de Marte.
        const double mars_theta_updated = mars_theta + MARS_INIT_V_THETA * DELTA_T;

        //      Posição radial da Sonda.
        const double ship_r_updated = ship_r + ship_v * DELTA_T;

        //      Correção da velocidade radial da Sonda.
        const double ship_v_updated = ship_v +
            ((ship_r * ship_v_theta * ship_v_theta) - (G_CONST * SUN_M / (ship_r * ship_r)) +
                (G_CONST * MARS_R * (ship_r - MARS_R * cos(ship_theta - mars_theta)) / (div * sqrt(div)))) * DELTA_T;

        //      Posição angular da Sonda.
        const double ship_theta_updated = ship_theta + ship_v_theta * DELTA_T;

        //      Correção da velocidade angular da Sonda.
        const double ship_v_theta_updated = ship_v_theta + ((G_CONST * MARS_M * MARS_R * sin(ship_theta - mars_theta) /
            (ship_r * div * sqrt(div))) - 2 * ship_v * ship_v_theta / ship_r) * DELTA_T;

        mars_theta = mars_theta_updated;
        ship_r = ship_r_updated;
        ship_v = ship_v_updated;
        ship_theta = ship_theta_updated;
        ship_v_theta = ship_v_theta_updated;
    }

    const double v_out = module_v(ship_v, ship_theta, ship_r, ship_v_theta);

    *rel_velocity = v_out - v_in;
    *deflection = comp_deflection(theta_init_value, ship_r, ship_v, ship_theta, ship_v_theta, mars_theta);

    fclose(fo);
}

int main(void) {
    double b_values[NUMBER_OF_TESTS];           //  [m]         -   Parâmetros de impacto usados no teste.
    double theta_init_ship[NUMBER_OF_TESTS];    //  [rad]       -   Ângulo inicial da nave.
    double delta_velocity[NUMBER_OF_TESTS];     //  [m/s]       -   Variação de velocidade entre a entreda e a saída da nave.
    double deflection_angle[NUMBER_OF_TESTS];   //  [rad]       -   Ângulo de deflexão.

    const double b_step = (B_MAX_VALUE - B_MIN_VALUE) / (NUMBER_OF_TESTS - 1);
    b_values[0] = B_MIN_VALUE;
    printf("b values: %.3e", b_values[0]);
    for (int i = 1; i < NUMBER_OF_TESTS; i++) {
        b_values[i] = b_values[i - 1] + b_step;
        printf("\t%.3e", b_values[i]);
    }

    printf("\n");
    printf("theta ship initial values: ");
    for (int i = 0; i < NUMBER_OF_TESTS; i++) {
        theta_init_ship[i] = MARS_INIT_THETA * (M_PI / 180.0) + asin(b_values[i] / SHIP_R);
        printf("\t%.6e", theta_init_ship[i]);
    }

    printf("\n");
    for (int i = 0; i < NUMBER_OF_TESTS; i++)
        simulate(i, &delta_velocity[i], &deflection_angle[i], theta_init_ship[i]);

    FILE *fo = fopen("../data/global.dat", "w");
    for (int i = 0; i < NUMBER_OF_TESTS; i++)
        fprintf(fo, "%e %e %e\n", b_values[i], delta_velocity[i], deflection_angle[i]);

    fclose(fo);

    return 0;
}