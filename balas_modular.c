#include "balas_modules.h"

int main(int argc, char *argv[]){
    // Read problem data
    read_problem_data(argv[1]);
    float max_time_seconds = atof(argv[2]);
    clock_t end_time;
    clock_t start_time = clock();
    clock_t stopping_time = start_time + max_time_seconds * CLOCKS_PER_SEC;

    initialize_infeasible_rows(constraint_matrix, row_starts, variable_indices, right_hand_sides, number_of_constraints);

    // if there are no infeasible rows initially, then the solution is just x1 = x2 = ... = xN = 0
    if (head_row == NULL){
        printf("Optimal solution is x1 = x2 = ... = xN = 0.\n");
        return 0;
    }

    initialize_free_variables(objective_coefficients, number_of_variables);
    while (!algorithm_done && clock() < stopping_time){
        execute_iteration();
    }

    if (best_objective_value == DBL_MAX){
        end_time = clock();
        printf("Problem is infeasible.\n");
        printf("Total execution time: %.6lf s.\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);
        return -1;
    }

    end_time = clock();
    printf("Total execution time: %.6lf s.\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);
    print_best_found_solution();

    free_all_memory();
    return 0;
}