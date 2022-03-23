#include "balas_modules.h"

int main(){
    clock_t start_time = clock();
    clock_t end_time;

    // TODO implement better I/O

    initialize_infeasible_rows(constraint_matrix, row_starts, variable_indices, right_hand_sides, number_of_constraints);

    // if there are no infeasible rows initially, then the solution is just x1 = x2 = ... = xN = 0
    if (head_row == NULL){
        printf("Optimal solution is x1 = x2 = ... = xN = 0.\n");
        end_time = clock();
        printf("Total execution time: %.6lf s.\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);
        return 0;
    }

    initialize_free_variables(objective_coefficents, number_of_variables);
    while (!algorithm_done){
        execute_iteration();
    }

    if (best_objective_value == DBL_MAX){
        printf("Problem is infeasible.\n");
        end_time = clock();
        printf("Total execution time: %.6lf s.\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);
        return -1;
    }

    // TODO implement solution printing

    end_time = clock();
    printf("Total execution time: %.6lf s.\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);
    return 0;
}