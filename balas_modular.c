#include "balas_modules.h"
#include "problem_data.h"

int main(){
    initialize_infeasible_rows(constraint_matrix, row_starts, variable_indices, right_hand_sides, number_of_constraints);

    // if there are no infeasible rows initially, then the solution is just x1 = x2 = ... = xN = 0
    if (head_row == NULL){
        printf("Optimal solution is x1 = x2 = ... = xN = 0.\n");
        return 0;
    }

    initialize_free_variables(objective_coefficents, number_of_variables);
    while(!algorithm_done){
        // for each free variable, find infeasible rows containing that variable
    }

    return 0;
}