#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <time.h>
#include "problem_data.h"

bool algorithm_done = false;
double current_objective_value = 0;
double best_objective_value = DBL_MAX;
bool *current_solution = incumbent_solution;

struct infeasible_row{
    int row_index;
    double *coefficients;
    int *variable_indices;
    int number_of_nonzeros;
    struct infeasible_row *next;
};

struct infeasible_row *head_row = NULL;

struct free_variable{
    int index;
    double objective_coefficient;
    struct free_variable *next;
};

struct free_variable *head_free_variable = NULL;

struct infeasibility_reducing_variable{
    int index;
    struct infeasibility_reducing_variable *next;
};

struct infeasibility_reducing_variable *head_infeasibility_reducing_variable = NULL;

struct fixed_variable{
    int index;
    bool fixed_to_one;
    double objective_coefficient;
    struct fixed_variable *next;
};

struct fixed_variable *head_fixed_variable = NULL;

void insert_row_front(int row_index, double *coefficients, int *variable_indices, int number_of_nonzeros){
    /*
    Inserts an infeasible row at the front of the list of infeasible rows
    Parameters: // TODO
    */
    struct infeasible_row *row = (struct infeasible_row*) malloc(sizeof(struct infeasible_row));
    row->row_index = row_index;
    row->coefficients = coefficients;
    row->variable_indices = variable_indices;
    row->number_of_nonzeros = number_of_nonzeros;
    row->next = head_row;
    head_row = row;
}

void delete_infeasible_row(int row_index){
    /*
    Deletes an infeasible row, specified by the row start.
    This method has no return value, and doesn't modify the list if
    the requested row is not found or if the list is already empty.
    */
    struct infeasible_row *current_row = head_row;
    struct infeasible_row *previous_row = NULL;

    // if list is empty
    if (head_row == NULL){
        return;
    }

    // go through list starting at first row to find the requested row
    while (current_row->row_index != row_index){
        // if we're at the last row
        if (current_row->next == NULL){
            return;
        }
        else{
            previous_row = current_row;
            current_row = current_row->next;
        }
    }

    if (current_row == head_row){
        head_row = head_row->next;
    }
    else{
        previous_row->next = current_row->next;
    }

    free(current_row);
}

void delete_free_variable(int index){
    /*
    Deletes a free variable, specified by its index.
    This method has no return value, and doesn't modify the list if
    the requested variable is not found or if the list is already empty.
    */
    struct free_variable *current_variable = head_free_variable;
    struct free_variable *previous_variable = NULL;

    // if list is empty
    if (head_free_variable == NULL){
        return;
    }

    // go through list starting at first variable to find the requested variable
    while (current_variable->index != index){
        // if we're at the last variable
        if (current_variable->next == NULL){
            return;
        }
        else{
            previous_variable = current_variable;
            current_variable = current_variable->next;
        }
    }

    if (current_variable == head_free_variable){
        head_free_variable = head_free_variable->next;
    }
    else{
        previous_variable->next = current_variable->next;
    }

    free(current_variable);
}

void insert_free_variable_front(int index, double objective_coefficient){
    struct free_variable *variable = (struct free_variable*) malloc(sizeof(struct free_variable));
    variable->index = index;
    variable->objective_coefficient = objective_coefficient;
    variable->next = head_free_variable;
    head_free_variable = variable;
}

void insert_infeasibility_reducing_variable_front(int index){
    struct infeasibility_reducing_variable *variable = (struct infeasibility_reducing_variable*) malloc(sizeof(struct infeasibility_reducing_variable));
    variable->index = index;
    variable->next = head_infeasibility_reducing_variable;
    head_infeasibility_reducing_variable = variable;
}

void insert_fixed_variable_front(int index, double objective_coefficient){
    struct fixed_variable *variable = (struct fixed_variable*) malloc(sizeof(struct fixed_variable));
    variable->index = index;
    variable->fixed_to_one = true;
    variable->objective_coefficient = objective_coefficient;
    variable->next = head_fixed_variable;
    head_fixed_variable = variable;
}

double calculate_objective_value(){
    double objective_value = 0;
    struct fixed_variable *current_fixed_variable = head_fixed_variable;
    while(current_fixed_variable != NULL){
        if (current_fixed_variable->fixed_to_one){
            objective_value += current_fixed_variable->objective_coefficient;
        }
        current_fixed_variable = current_fixed_variable->next;
    }
    return objective_value;
}

bool row_contains_free_variable(struct infeasible_row *row, struct free_variable *variable){
    for (int i = 0; i < row->number_of_nonzeros; i++){
        if (row->variable_indices[i] == variable->index){
            return true;
        }
    }
    return false;
}

bool row_contains_infeasibility_reducing_variable(struct infeasible_row *row, struct infeasibility_reducing_variable *variable){
    for (int i = 0; i < row->number_of_nonzeros; i++){
        if (row->variable_indices[i] == variable->index){
            return true;
        }
    }
    return false;
}

bool row_is_infeasible(int row_index){
    struct infeasible_row *current_row = head_row;
    while (current_row != NULL){
        if (current_row->row_index == row_index){
            return true;
        }
        current_row = current_row->next;
    }
    return false;
}

double calculate_infeasibility_reduction(struct infeasible_row *row, struct infeasibility_reducing_variable *variable){
    for (int i = 0; i < row->number_of_nonzeros; i++){
        if ((row->variable_indices[i] == variable->index) && (row->coefficients[i] < 0)){
            return row->coefficients[i];
        }
    }
}

double calculate_minimum_infeasibility(struct infeasible_row *row){
    double minimum_infeasibility = left_hand_sides[row->row_index] - right_hand_sides[row->row_index];
    struct infeasibility_reducing_variable *current_infeasibility_reducing_variable = head_infeasibility_reducing_variable;
    while (current_infeasibility_reducing_variable != NULL){
        if (row_contains_infeasibility_reducing_variable(row, current_infeasibility_reducing_variable)){
            double infeasibility_reduction = calculate_infeasibility_reduction(row, current_infeasibility_reducing_variable);
            minimum_infeasibility += infeasibility_reduction; // infeasibility_reduction is a negative value
        }
        current_infeasibility_reducing_variable = current_infeasibility_reducing_variable->next;
    }
    return minimum_infeasibility;
}

void delete_all_infeasibility_reducing_variables(){
    struct infeasibility_reducing_variable *current_infeasibility_reducing_variable = head_infeasibility_reducing_variable;
    struct infeasibility_reducing_variable *next_infeasibility_reducing_variable = NULL;
    while (current_infeasibility_reducing_variable != NULL){
        next_infeasibility_reducing_variable = current_infeasibility_reducing_variable->next;
        free(current_infeasibility_reducing_variable);
        current_infeasibility_reducing_variable = next_infeasibility_reducing_variable;
    }
    head_infeasibility_reducing_variable = NULL;
}

void initialize_infeasible_rows(double *constraint_matrix, int *row_starts, int *variable_indices, double *right_hand_sides, int number_of_constraints){
    /*
    At the start of the algorithm, all variable are free and evaluated at zero
    Since we assume the constraint matrix is of the form Ax <= b, the i-th row is initially infeasible
    if and only if b_i < 0
    */
    for (int i = 0; i < number_of_constraints; i++){
        if (*(right_hand_sides + i) < 0){
            int row_start = *(row_starts + i);
            int number_of_nonzeros = *(row_starts + i + 1) - row_start;
            insert_row_front(row_start, constraint_matrix + row_start, variable_indices + row_start, number_of_nonzeros);
        }
    }
    return;
}

void initialize_free_variables(double *objective_coefficients, int number_of_variables){
    /*
    Initially, every variable is free
    */
    for (int i = 0; i < number_of_variables; i++){
        insert_free_variable_front(i, *(objective_coefficients + i));
    }
}

void update_infeasible_rows(int variable_index, int direction){
    /*
    When a variable is newly fixed to 0 or 1, we recalculate the left-hand side for all rows accordingly
    and update the list of infeasible rows, adding or removing rows as needed. The input 'direction' is 0
    or 1, indicating which value the variable has been fixed to.
    */
    for (int row_index = 0; row_index < number_of_constraints; row_index++){
        for (int i = row_starts[row_index]; i < row_starts[row_index + 1]; i++){
            if (variable_indices[i] == variable_index){
                if (direction){
                    left_hand_sides[row_index] += constraint_matrix[i];
                }
                else{
                    left_hand_sides[row_index] -= constraint_matrix[i];
                }
                // update lists
                if (row_is_infeasible(row_index)){
                    // If the row is now feasible, delete it from the list of infeasible rows
                    if (left_hand_sides[row_index] <= right_hand_sides[row_index]){
                        delete_infeasible_row(row_index);
                    }
                }
                else{
                    // If the row is now infeasible, add it to the list of infeasible rows
                    if (left_hand_sides[row_index] > right_hand_sides[row_index]){
                        insert_row_front(row_index, constraint_matrix + row_starts[row_index], variable_indices + row_starts[row_index], 
                            row_starts[row_index + 1] - row_starts[row_index]
                        );
                    }
                }

                // move to the next constraint
                break;
            }
        }
    }
}

void backtrack(){
    /*
    Backtracking means:
    - if the last fixed variable is fixed to 1, set it to 0
    - if the last fixed variable is fixed to 0, change it to a free variable and consider the second-to-last fixed variable
        - if that variable is fixed to 1, set it to 0
        - if it's fixed to 0, change it to a free variable and backtrack one further
    - if we backtrack all the way to unfixing the first fixed variable, the algorithm is done running
    */

    // if no fixed variables remain, we are done with Balas' algorithm
    if (head_fixed_variable == NULL){
        algorithm_done = true;
        return;
    }

    // otherwise, proceed with backtracking
    struct fixed_variable *current_fixed_variable = head_fixed_variable;
    if (current_fixed_variable->fixed_to_one){
        // fix to zero and modify the current objective value
        current_fixed_variable->fixed_to_one = false;
        update_infeasible_rows(current_fixed_variable->index, 0);
        current_objective_value -= current_fixed_variable->objective_coefficient;
        return;
    }
    else{
        // if fixed to zero, change head fixed variable, unfix current fixed variable and add to list of free variables
        head_fixed_variable = current_fixed_variable->next;
        insert_free_variable_front(current_fixed_variable->index, current_fixed_variable->objective_coefficient);
        free(current_fixed_variable);
        backtrack();
    }
}

double calculate_total_infeasibility_reduction(struct free_variable *variable){
    /*
    Calculates the amount by which setting a free variable to 1 will reduce the total infeasibility
    of infeasible rows
    */
    struct infeasible_row *current_row = head_row;
    double total_infeasibility_reduction = 0;
    while (current_row != NULL){
        for (int i = 0; i < current_row->number_of_nonzeros; i++){
            if ((current_row->variable_indices[i] == variable->index) && (current_row->coefficients[i] < 0)){
                total_infeasibility_reduction += current_row->coefficients[i];
            }
        }
        current_row = current_row->next;
    }
    return total_infeasibility_reduction;
}

void update_incumbent_solution(){
    /*
    Stores the current solution in the incumbent slot.
    */
    for (int i = 0; i < number_of_variables; i++){
        incumbent_solution[i] = false;
    }
    struct fixed_variable *current_fixed_variable = head_fixed_variable;
    while (current_fixed_variable != NULL){
        if (current_fixed_variable->fixed_to_one){
            incumbent_solution[current_fixed_variable->index] = true;
        }
        current_fixed_variable = current_fixed_variable->next;
    }
}

void execute_iteration(){
    /*
    If there are no infeasible rows left, we can prune this node 
    */
    if (head_row == NULL){
        if (current_objective_value < best_objective_value){
            update_incumbent_solution();
        }
        return;
    }

    /*
    For each free variable, check if setting it to 1 would reduce infeasibility and lead to a better objective than the incumbent
    The set of infeasibility reducing variables is the set T in the book with the PASCAL implementation
    */
    struct free_variable *current_free_variable = head_free_variable;
    while (current_free_variable != NULL){
        if (current_objective_value + current_free_variable->objective_coefficient < best_objective_value){
            struct infeasible_row *current_row = head_row;
            while (current_row != NULL){
                if (row_contains_free_variable(current_row, current_free_variable)){
                    /*
                    If a free variable reduces infeasibility in at least one row, that's enough to include it 
                    in the list of infeasibility reducing variables and we don't need to consider what other
                    rows it may appear in.
                    */
                    insert_infeasibility_reducing_variable_front(current_free_variable->index);
                    break;
                }
                current_row = current_row->next;
            }
        }
        current_free_variable = current_free_variable->next;
    }

    if (head_infeasibility_reducing_variable == NULL){
        // If there are no infeasibility reducing variables (T is empty), backtrack
        backtrack();
        return;
    }

    /*
    Proceed with infeasibility test (if we set all infeasibility reducing variables to 1, will the
    solution still be infeasible?)
    We check every infeasible row
    */
    struct infeasible_row *current_row = head_row;
    while (current_row != NULL){
        double minimum_infeasibility = calculate_minimum_infeasibility(current_row);
        if (minimum_infeasibility > 0){
            /*
            If there exists a row that cannot be made feasible by setting all infeasibility reducing variables 
            to 1, then there is no feasible completion of the current partial solution, so we backtrack
            */
            backtrack();
            return;
        }
        current_row = current_row->next;
    }

    // We no longer need the list of infeasibility reducing variables
    delete_all_infeasibility_reducing_variables();
    
    /*
    There is a feasible continuation of the partial solution. We decide which variable to branch on using Balas' test.
    By convention, infeasibility reductions are negative if setting a variable to 1 reduces the infeasibility. 
    */
    double largest_infeasibility_reduction = 0;
    int branching_variable_index;
    double branching_variable_objective_coefficient;
    current_free_variable = head_free_variable;
    while (current_free_variable != NULL){
        double infeasibility_reduction = calculate_total_infeasibility_reduction(current_free_variable);
        if (infeasibility_reduction < largest_infeasibility_reduction){
            largest_infeasibility_reduction = infeasibility_reduction;
            branching_variable_index = current_free_variable->index;
            branching_variable_objective_coefficient = current_free_variable->objective_coefficient;
        }
        current_free_variable = current_free_variable->next;
    }
    // If the largest infeasibility reduction is 0, then we cannot reduce infeasibility by setting a variable to 1
    if (largest_infeasibility_reduction == 0){
        backtrack();
        return;
    }
    // Branch on the variable with the largest infeasibility reduction
    insert_fixed_variable_front(branching_variable_index, branching_variable_objective_coefficient);
    delete_free_variable(branching_variable_index);
    update_infeasible_rows(branching_variable_index, 1);
    current_objective_value += branching_variable_objective_coefficient;
    if (current_objective_value < best_objective_value){
        best_objective_value = current_objective_value;
        update_incumbent_solution();
    }
}

void print_solution(){
    // TODO make this more sensible
    printf("Optimal solution:\n[");
    for (int i = 0; i < number_of_variables - 1; i++){
        if (incumbent_solution[i]){
            printf("1");
        }
        else{
            printf("0");
        }
        printf(", ");
    }
    if (incumbent_solution[number_of_variables - 1]){
        printf("1");
    }
    else{
        printf("0");
    }
    printf("]\n");
    printf("Optimal objective value: %.5lf\n", best_objective_value);
}