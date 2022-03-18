#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

struct infeasible_row{
    int index;
    double infeasibility;
    double *coefficients;
    int *variable_indices;
    int number_of_nonzeros;
    struct infeasible_row *next;
};

struct infeasible_row *head_row = NULL;
struct infeasible_row *current_row = NULL;

struct free_variable{
    int index;
    double objective_coefficient;
    struct free_variable *next;
};

struct free_variable *head_free_variable = NULL;

void insert_row_front(int index, double infeasibility, double *coefficients, int *variable_indices, int number_of_nonzeros){
    /*
    Inserts an infeasible row at the front of the list of infeasible rows
    Parameters:
        index: index of infeasible row in the constraint matrix
        infeasibility: amount by which row is currently infeasible
    */
    struct infeasible_row *row = (struct infeasible_row*) malloc(sizeof(struct infeasible_row));
    row->index = index;
    row->infeasibility = infeasibility;
    row->coefficients = coefficients;
    row->variable_indices = variable_indices;
    row->number_of_nonzeros = number_of_nonzeros;
    row->next = head_row;
    head_row = row;
}

void delete_infeasible_row(int index){
    /*
    Deletes an infeasible row, specified by a row index.
    This method has no return value, and doesn't modify the list if
    the requested row is not found or if the list is already empty.
    */
    struct infeasible_row *current_row = head_row;
    struct infeasible_row *previous = NULL;

    // if list is empty
    if (head_row == NULL){
        return;
    }

    // go through list starting at first row to find the requested row
    while (current_row->index != index){
        // if we're at the last row
        if (current_row->next == NULL){
            return;
        }
        else{
            previous = current_row;
            current_row = current_row->next;
        }
    }

    if (current_row == head_row){
        head_row = head_row->next;
    }
    else{
        previous->next = current_row->next;
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

double calculate_objective_value(double *coefficients, bool *solution, int number_of_variables){
    double objective_value = 0;
    for (int i = 0; i < number_of_variables; i++){
        if (*(solution + i)){
            objective_value += coefficients[i];
        }
    }
    return objective_value;
}

bool row_contains_variable(struct infeasible_row *row, struct free_variable *free_variable){
    for (int i = 0; i < row->number_of_nonzeros; i++){
        if (row->variable_indices[i] == free_variable->index){
            return true;
        }
    }
    return false;
}

double calculate_infeasibility_reduction(struct infeasible_row *row, struct free_variable *variable){ // TODO
    return 0;
}

double calculate_minimum_infeasibility(struct infeasible_row *row){
    double minimum_infeasibility = row->infeasibility;
    struct free_variable *current_free_variable = head_free_variable;
    while (current_free_variable != NULL){
        if (row_contains_variable(row, current_free_variable)){
            double infeasibility_reduction = calculate_infeasibility_reduction(row, current_free_variable);
            minimum_infeasibility -= infeasibility_reduction;
        }
    }
    return minimum_infeasibility;
}

void initialize_infeasible_rows(){ //TODO
    return;
}

void initialize_free_variables(){ //TODO
    return;
}