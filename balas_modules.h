#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

bool algorithm_done = false;

struct infeasible_row{
    int row_start;
    double infeasibility;
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

void insert_row_front(int row_start, double infeasibility, double *coefficients, int *variable_indices, int number_of_nonzeros){
    /*
    Inserts an infeasible row at the front of the list of infeasible rows
    Parameters: // TODO
    */
    struct infeasible_row *row = (struct infeasible_row*) malloc(sizeof(struct infeasible_row));
    row->row_start = row_start;
    row->infeasibility = infeasibility;
    row->coefficients = coefficients;
    row->variable_indices = variable_indices;
    row->number_of_nonzeros = number_of_nonzeros;
    row->next = head_row;
    head_row = row;
}

void delete_infeasible_row(int row_start){
    /*
    Deletes an infeasible row, specified by a the row start.
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
    while (current_row->row_start != row_start){
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

double calculate_row_infeasibility(bool *solution, double *constraint_matrix, int *constraint_columns, int row_start, int row_length, double row_rhs){
    /*
    Given a constraint matrix in CSR format, Ax <= b, for a binary IP, this function calculates the amount by which a constraint is violated
    Parameters:
        solution: current solution vector [0,1,0,1,1,...]
        solution_length: number of variables (N)
        constraint_matrix: array of CSR matrix coefficients
        constraint_columns: array of CSR matrix variable indices
        row_start: index of the start of the constraint in question
        row_length: number of nonzeroes in the constraint in question
        row_rhs: right-hand side of the row (b_i)
    Returns:
        infeasibility: amount by which the constraint is violated, i.e. sum(i, a_i * x_i) - b_i
            if infeasibility is positive, it means the solution violates the constraint
            if it's zero or negative, solution satisfies the constraint
    */
   double infeasibility = -row_rhs;
   for (int i = 0; i < row_length; i++){
       int variable_index = *(constraint_columns + row_start + i);
       if (*(solution + variable_index)){
           infeasibility += (*(constraint_matrix + row_start + i));
       }
   }
   return infeasibility;
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

double calculate_infeasibility_reduction(struct infeasible_row *row, struct infeasibility_reducing_variable *variable){
    for (int i = 0; i < row->number_of_nonzeros; i++){
        if (row->variable_indices[i] == variable->index){
            return row->coefficients[i];
        }
    }
}

double calculate_minimum_infeasibility(struct infeasible_row *row){
    double minimum_infeasibility = row->infeasibility;
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
            insert_row_front(row_start, -*(right_hand_sides + i), constraint_matrix + row_start, variable_indices + row_start, number_of_nonzeros);
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
        // fix to zero
        current_fixed_variable->fixed_to_one = false;
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