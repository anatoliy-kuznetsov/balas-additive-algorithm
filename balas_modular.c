#include "balas_modules.h"
#include "problem_data.h"

int main(){
    initialize_infeasible_rows(constraint_matrix, row_starts, variable_indices, right_hand_sides, number_of_constraints);
    return 0;
}