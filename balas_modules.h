#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

struct infeasible_row{
    int index;
    double infeasibility;
    struct infeasible_row *next;
};

struct infeasible_row *head_row = NULL;
struct infeasible_row *current_row = NULL;

void insert_row_front(int index, double infeasibility){
    /*
    Inserts an infeasible row at the front of the list of infeasible rows
    Parameters:
        index: index of infeasible row in the constraint matrix
        infeasibility: amount by which row is infeasible
    */
    struct infeasible_row *row = (struct infeasible_row*) malloc(sizeof(struct infeasible_row));
    row->index = index;
    row->infeasibility = infeasibility;
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
}

double calculate_row_infeasibility(bool *solution, int solution_length, double *constraint_matrix, int *constraint_columns, int row_start, int row_length, double row_rhs){
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