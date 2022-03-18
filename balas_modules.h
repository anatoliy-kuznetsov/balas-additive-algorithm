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