#include "balas_modules.h"

int main(int argc, char *argv[]){
    /*
    Command-line arguments:
        - argv[1] is the name of the input file
        - argv[2] is the time limit in seconds (double)
    Example:
    ./balas_modular balas-example1.balasin 0.5
    */
    read_problem_data(argv[1]);
    double max_time_seconds = atof(argv[2]);
    clock_t end_time;
    clock_t start_time = clock();
    clock_t stopping_time = start_time + max_time_seconds * CLOCKS_PER_SEC;

    initialize_infeasible_rows(constraint_matrix, row_starts, variable_indices, right_hand_sides, number_of_constraints);

    // if there are no infeasible rows initially, then the solution is just x1 = x2 = ... = xN = 0
    if (head_row == NULL){
        printf("Optimal solution is x1 = x2 = ... = xN = 0.\n");
        free_all_memory();
        return 0;
    }

    initialize_free_variables(objective_coefficients, number_of_variables);
    while (!optimal_solution_found && clock() < stopping_time){
        execute_iteration();
        iteration_counter++;
    }

    if (best_objective_value == DBL_MAX){
        end_time = clock();
        if (end_time >= stopping_time){
            printf("Maximum time exceeded. No feasible solution found.\n");
        }
        else{
            printf("Problem is infeasible.\n");
        }
        printf("Total execution time: %.6lf s.\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);
        printf("Iterations: %lld, time per iteration: %.9lf s.\n", iteration_counter, (double)(end_time - start_time)/(CLOCKS_PER_SEC * iteration_counter));
        free_all_memory();
        return -1;
    }

    end_time = clock();
    printf("Total execution time: %.6lf s.\n", (double)(end_time - start_time)/CLOCKS_PER_SEC);
    printf("Iterations: %lld, time per iteration: %.9lf s.\n", iteration_counter, (double)(end_time - start_time)/(CLOCKS_PER_SEC * iteration_counter));
    print_best_found_solution();

    free_all_memory();
    return 0;
}