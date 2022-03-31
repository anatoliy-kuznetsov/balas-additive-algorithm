#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <float.h>
#include <time.h>

bool optimal_solution_found = false;
double current_objective_value = 0;
double best_objective_value = DBL_MAX;
bool *incumbent_solution;
bool *current_solution;

int number_of_variables;
int number_of_constraints;
int number_of_nonzeros;
int number_of_objective_nonzeros;
double *constraint_matrix;
int *row_starts;
double *objective_coefficients;
int *objective_indices;
int *variable_indices;
double *right_hand_sides;
double *left_hand_sides;

struct infeasible_row{
    int row_index;
    double *coefficients;
    int *variable_indices;
    int number_of_nonzeros;
    double minimum_infeasibility;
    struct infeasible_row *next;
};

struct infeasible_row *head_row = NULL;

struct free_variable{
    int index;
    double objective_coefficient;
    double infeasibility_reduction;
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

struct row_with_negative_coefficient{
    int row_index;
    double coefficient;
    struct row_with_negative_coefficient *next;
};

struct row_with_nonzero_coefficient{
    int row_index;
    double coefficient;
    struct row_with_nonzero_coefficient *next;
};

struct row_with_negative_coefficient **rows_with_negative_coefficients;
struct row_with_nonzero_coefficient **rows_with_nonzero_coefficients;

void insert_row_with_nonzero_coefficient_end(int variable_index, int row_index, double coefficient){
    /*
    Given a variable and the index of a row that contains that variable, inserts a node at the end of the 
    corresponding list of rows that contain that variable
    */
    // If list is empty (head pointer is null) then we put the node at the start of the list
    if (rows_with_nonzero_coefficients[variable_index] == NULL){
        rows_with_nonzero_coefficients[variable_index] = (struct row_with_nonzero_coefficient*) malloc(sizeof(struct row_with_nonzero_coefficient));
        rows_with_nonzero_coefficients[variable_index]->row_index = row_index;
        rows_with_nonzero_coefficients[variable_index]->coefficient = coefficient;
        rows_with_nonzero_coefficients[variable_index]->next = NULL;
        return;
    }
    // Otherwise, traverse the list until we get to the end
    struct row_with_nonzero_coefficient *current_row_with_nonzero_coefficient = rows_with_nonzero_coefficients[variable_index];
    while (current_row_with_nonzero_coefficient->next != NULL){
        current_row_with_nonzero_coefficient = current_row_with_nonzero_coefficient->next;
    }
    struct row_with_nonzero_coefficient *new_row_with_nonzero_coefficient = (struct row_with_nonzero_coefficient*) malloc(sizeof(struct row_with_nonzero_coefficient));
    new_row_with_nonzero_coefficient->row_index = row_index;
    new_row_with_nonzero_coefficient->coefficient = coefficient;
    new_row_with_nonzero_coefficient->next = NULL;
    current_row_with_nonzero_coefficient->next = new_row_with_nonzero_coefficient;
}

void insert_row_with_negative_coefficient_end(int variable_index, int row_index, double coefficient){
    /*
    Given a variable and the index of a row that contains that variable with a negative coefficient,
    inserts a node at the end of the corresponding list of rows that contain that variable
    */
    // If list is empty (head pointer is null) then we put the node at the start of the list
    if (rows_with_negative_coefficients[variable_index] == NULL){
        rows_with_negative_coefficients[variable_index] = (struct row_with_negative_coefficient*) malloc(sizeof(struct row_with_negative_coefficient));
        rows_with_negative_coefficients[variable_index]->row_index = row_index;
        rows_with_negative_coefficients[variable_index]->coefficient = coefficient;
        rows_with_negative_coefficients[variable_index]->next = NULL;
        return;
    }
    // Otherwise, traverse the list until we get to the end
    struct row_with_negative_coefficient *current_row_with_negative_coefficient = rows_with_negative_coefficients[variable_index];
    while (current_row_with_negative_coefficient->next != NULL){
        current_row_with_negative_coefficient = current_row_with_negative_coefficient->next;
    }
    struct row_with_negative_coefficient *new_row_with_negative_coefficient = (struct row_with_negative_coefficient*) malloc(sizeof(struct row_with_negative_coefficient));
    new_row_with_negative_coefficient->row_index = row_index;
    new_row_with_negative_coefficient->coefficient = coefficient;
    new_row_with_negative_coefficient->next = NULL;
    current_row_with_negative_coefficient->next = new_row_with_negative_coefficient;
}

void insert_row_front(int row_index, double *coefficients, int *variable_indices, int number_of_nonzeros){
    /*
    Inserts an infeasible row at the front of the list of infeasible rows
    Parameters:
        row_index: which row this is in the constraint matrix
        coefficients: CSR representation of coefficients in this row
        variable_indices: array of nonzero column indices
        number_of_nonzeros: number of nonzero column indices
    */
    struct infeasible_row *row = (struct infeasible_row*) malloc(sizeof(struct infeasible_row));
    row->row_index = row_index;
    row->coefficients = coefficients;
    row->variable_indices = variable_indices;
    row->number_of_nonzeros = number_of_nonzeros;
    row->minimum_infeasibility = left_hand_sides[row_index] - right_hand_sides[row_index];
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
    variable->infeasibility_reduction = 0;
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

bool row_contains_infeasibility_reducing_variable(struct infeasible_row *row, struct infeasibility_reducing_variable *variable){
    for (int i = 0; i < row->number_of_nonzeros; i++){
        if (row->variable_indices[i] == variable->index){
            return true;
        }
    }
    return false;
}

bool row_in_infeasible_list(int row_index){
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
    return 0;
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

void delete_all_free_variables(){
    struct free_variable *current_free_variable = head_free_variable;
    struct free_variable *next_free_variable = NULL;
    while (current_free_variable != NULL){
        next_free_variable = current_free_variable->next;
        free(current_free_variable);
        current_free_variable = next_free_variable;
    }
    head_free_variable = NULL;
}

void delete_all_infeasible_rows(){
    struct infeasible_row *current_row = head_row;
    struct infeasible_row *next_row = NULL;
    while (current_row != NULL){
        next_row = current_row->next;
        free(current_row);
        current_row = next_row;
    }
    head_row = NULL;
}

void delete_all_row_lists(){
    for (int i = 0; i < number_of_variables; i++){
        struct row_with_nonzero_coefficient *current_row_with_nonzero_coefficient = rows_with_nonzero_coefficients[i];
        struct row_with_negative_coefficient *current_row_with_negative_coefficient = rows_with_negative_coefficients[i];
        struct row_with_nonzero_coefficient *next_row_with_nonzero_coefficient = NULL;
        struct row_with_negative_coefficient *next_row_with_negative_coefficient = NULL;
        while (current_row_with_nonzero_coefficient != NULL){
            next_row_with_nonzero_coefficient = current_row_with_nonzero_coefficient->next;
            free(current_row_with_nonzero_coefficient);
            current_row_with_nonzero_coefficient = next_row_with_nonzero_coefficient;
        }
        while (current_row_with_negative_coefficient != NULL){
            next_row_with_negative_coefficient = current_row_with_negative_coefficient->next;
            free(current_row_with_negative_coefficient);
            current_row_with_negative_coefficient = next_row_with_negative_coefficient;
        }
        rows_with_nonzero_coefficients[i] = NULL;
        rows_with_negative_coefficients[i] = NULL;
    }
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
            insert_row_front(i, constraint_matrix + row_start, variable_indices + row_start, number_of_nonzeros);
        }
    }
    return;
}

int get_objective_index(int variable_index){
    /*
    Returns the index in the sparse objective function of a variable with original index variable_index
    If the variable does not appear in the objective function, returns -1.
    */
    for (int i = 0; i < number_of_objective_nonzeros; i++){
        if (variable_index == objective_indices[i]){
            return i;
        }
    }
    return -1;
}

void initialize_free_variables(double *objective_coefficients, int number_of_variables){
    /*
    Initially, every variable is free
    */
    for (int i = 0; i < number_of_variables; i++){
        int objective_index = get_objective_index(i);
        int objective_coefficient;
        if (objective_index == -1){
            objective_coefficient = 0;
        }
        else{
            objective_coefficient = *(objective_coefficients + objective_index);
        }
        insert_free_variable_front(i, objective_coefficient);
    }
}

void update_infeasible_rows(int variable_index, int direction){
    /*
    When a variable is newly fixed to 0 or 1, we recalculate the left-hand side for all rows accordingly
    and update the list of infeasible rows, adding or removing rows as needed. The input 'direction' is 0
    or 1, indicating which value the variable has been fixed to.
    */
    struct row_with_nonzero_coefficient *current_row_with_nonzero_coefficient = rows_with_nonzero_coefficients[variable_index];
    while (current_row_with_nonzero_coefficient != NULL){
        int row_index = current_row_with_nonzero_coefficient->row_index;
        double previous_left_hand_side = left_hand_sides[row_index];
        double coefficient = current_row_with_nonzero_coefficient->coefficient;
        if (direction){
            left_hand_sides[row_index] += coefficient;
        }
        else{
            left_hand_sides[row_index] -= coefficient;
        }
        if (previous_left_hand_side > right_hand_sides[row_index] &&
            left_hand_sides[row_index] <= right_hand_sides[row_index]){
                // The row was made feasible
                delete_infeasible_row(row_index);
            }
        if (previous_left_hand_side <= right_hand_sides[row_index] &&
            left_hand_sides[row_index] > right_hand_sides[row_index]){
                // The row was made infeasible
                insert_row_front(row_index, constraint_matrix + row_starts[row_index], variable_indices + row_starts[row_index], 
                    row_starts[row_index + 1] - row_starts[row_index]
                );
            }
        current_row_with_nonzero_coefficient = current_row_with_nonzero_coefficient->next;
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
        optimal_solution_found = true;
        return;
    }

    // otherwise, proceed with backtracking
    struct fixed_variable *current_fixed_variable = head_fixed_variable;
    if (current_fixed_variable->fixed_to_one){
        // fix to zero, modify the current objective value, and update infeasible row list
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

void decrement_minimum_infeasibility(int row_index, double coefficient){
    /*
    Updates the value of minimum infeasibility for a row specified by its index.
    This method is called with the coefficient of an infeasibility reducing variable in that row
    */
    struct infeasible_row *current_row = head_row;
    while (current_row != NULL){
        if (current_row->row_index == row_index){
            current_row->minimum_infeasibility += coefficient;
            return;
        }
        current_row = current_row->next;
    }
}

bool is_free(int variable_index){
    /*
    Checks whether a variable, specified by its index, is currently free
    */
    struct free_variable *current_free_variable = head_free_variable;
    while (current_free_variable != NULL){
        if (current_free_variable->index == variable_index){
            return true;
        }
        current_free_variable = current_free_variable->next;
    }
    return false;
}

bool feasible_to_set_all_free_variables_to_one(){
    /*
    Checks whether setting all free variables to a value of 1 leads to a feasible solution.
    This method is called when there exists an infeasible row that can only be made feasible
    by setting all free variables to 1, implying there is only one possible feasible continuation of
    the partial solution. In this case, we need to check whether the rest of the rows will 
    remain feasible after such an operation.
    Checking this is rather slow, but most nodes in the tree probably don't meet this criterion
    */
    for (int row_index = 0; row_index < number_of_constraints; row_index++){
        double new_left_hand_side = left_hand_sides[row_index];
        for (int j = row_starts[row_index]; j < row_starts[row_index + 1]; j++){
            if (is_free(variable_indices[j])){
                new_left_hand_side += constraint_matrix[j];
            }
        }
        if (new_left_hand_side > right_hand_sides[row_index]){
            return false;
        }
    }
    return true;
}

void reset_minimum_infeasibilities(){
    /*
    At each iteration, we compute the minimum infeasibility of each infeasible row
    Between iterations, we set them all to [lhs] - [rhs]
    */
    struct infeasible_row *current_row = head_row;
    while (current_row != NULL){
        int row_index = current_row->row_index;
        current_row->minimum_infeasibility = left_hand_sides[row_index] - right_hand_sides[row_index];
        current_row = current_row->next;
    }
}

void reset_infeasibility_reductions(){
    /*
    At each iteration, we compute the infeasibility reduction of each variable
    Between iterations, we set them all to 0
    */
    struct free_variable *current_free_variable = head_free_variable;
    while (current_free_variable != NULL){
        current_free_variable->infeasibility_reduction = 0;
        current_free_variable = current_free_variable->next;
    }
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

void get_infeasibility_reducing_variables(double *largest_infeasibility_reduction, int *branching_variable_index, 
                                          double *branching_variable_objective_coefficient){
    /*
    Builds the set of infeasibility reducing variables and computes their infeasibility
    reduction, along with the minimum infeasibility of each infeasible row.
    Returns the largest infeasibility reduction of any variable, as well as the index and objective
    coefficient of the variable that reduces infeasibility the most
    */
    *largest_infeasibility_reduction = 0;
    struct free_variable *current_free_variable = head_free_variable;
    while (current_free_variable != NULL){
        bool variable_reduces_infeasibility = false;
        if (current_objective_value + current_free_variable->objective_coefficient < best_objective_value){
            struct row_with_negative_coefficient *current_row_with_negative_coefficient = rows_with_negative_coefficients[current_free_variable->index];
            while (current_row_with_negative_coefficient != NULL){
                int row_index = current_row_with_negative_coefficient->row_index;
                double coefficient = current_row_with_negative_coefficient->coefficient;
                /*
                If we haven't yet added the variable to the set of infeasibility reducing variables, we 
                do that and then set the corresponding flag to 1
                */
                if (!variable_reduces_infeasibility){
                    insert_infeasibility_reducing_variable_front(current_free_variable->index);
                    variable_reduces_infeasibility = true;
                }
                /*
                Calculate the total infeasibility reduction of the variable and add to the minimum infeasibility
                of each row it appears in with a negative coefficient
                */
                if (left_hand_sides[row_index] > right_hand_sides[row_index]){
                    current_free_variable->infeasibility_reduction += coefficient;
                    decrement_minimum_infeasibility(row_index, coefficient);
                }
                current_row_with_negative_coefficient = current_row_with_negative_coefficient->next;
            }
        }
        if (current_free_variable->infeasibility_reduction < *largest_infeasibility_reduction){
            *largest_infeasibility_reduction = current_free_variable->infeasibility_reduction;
            *branching_variable_index = current_free_variable->index;
            *branching_variable_objective_coefficient = current_free_variable->objective_coefficient;
        }
        current_free_variable = current_free_variable->next;
    }
}

double get_minimum_infeasibility_among_rows(){
    /*
    Carries out infeasibility test (if we set all infeasibility reducing variables to 1, will the
    solution still be infeasible?)
    We check every infeasible row
    */
    double minimum_infeasibility_among_rows = -DBL_MAX;
    struct infeasible_row *current_row = head_row;
    while (current_row != NULL){
        double minimum_infeasibility = current_row->minimum_infeasibility;
        if (minimum_infeasibility > 0){
            return minimum_infeasibility;
        }
        if (minimum_infeasibility > minimum_infeasibility_among_rows){
            minimum_infeasibility_among_rows = minimum_infeasibility;
        }
        current_row = current_row->next;
    }
    return minimum_infeasibility_among_rows;
}

void print_best_found_solution(){
    if (optimal_solution_found){
        printf("Optimal solution:\n");
    }
    else{
        printf("Maximum time exceeded. Best solution found:\n");
    }
    for (int i = 0; i < number_of_variables; i++){
        if (incumbent_solution[i]){
            printf("x%d:\t1\n", i + 1);
        }
        else{
            printf("x%d:\t0\n", i + 1);
        }
    }
    printf("Objective value: %.5lf\n", best_objective_value);
}

void execute_iteration(){
    /*
    Performs one iteration of Balas' algorithm
    */

    // If there are no infeasible rows left, we can prune this node 
    if (head_row == NULL){
        if (current_objective_value < best_objective_value){
            update_incumbent_solution();
            best_objective_value = current_objective_value;
        }
        backtrack();
        return;
    }
    
    double largest_infeasibility_reduction;
    int branching_variable_index;
    double branching_variable_objective_coefficient;
    get_infeasibility_reducing_variables(&largest_infeasibility_reduction, &branching_variable_index, &branching_variable_objective_coefficient);
    if (head_infeasibility_reducing_variable == NULL){
        // If there are no infeasibility reducing variables (the set T is empty), backtrack
        backtrack();
        return;
    }

    double minimum_infeasibility_among_rows = get_minimum_infeasibility_among_rows();
    if (minimum_infeasibility_among_rows > 0){
        /*
        If there exists a row that cannot be made feasible by setting all infeasibility reducing variables 
        to 1, then there is no feasible completion of the current partial solution, so we backtrack
        */
        backtrack();
        return;
    }

    /*
    If there is a feasible continuation and there is a row that can only be made feasible by setting all remaining variables to 1, 
    then there is only one possible feasible continuation and we know its objective value. A sufficient condition for this is for the 
    minimum infeasibility of any one row to be 0 (recall that a positive minimum infeasibility means a row cannot be made feasible, and
    a negative minimum infeasibility means that it can be made strictly feasible) and for the continuation to be feasible.
    */
    if (minimum_infeasibility_among_rows == 0){
        if (feasible_to_set_all_free_variables_to_one()){
            double only_possible_objective_value = current_objective_value;
            struct free_variable *current_free_variable = head_free_variable;
            while (current_free_variable != NULL){
                only_possible_objective_value += current_free_variable->objective_coefficient;
                current_free_variable = current_free_variable->next;
            }

            if (only_possible_objective_value < best_objective_value){
                best_objective_value = only_possible_objective_value;
                /*
                Store the solution by traversing the lists of free and fixed variables
                Similar to the method update_incumbent_solution() but with a small modification
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
                current_free_variable = head_free_variable;
                while (current_free_variable != NULL){
                    incumbent_solution[current_free_variable->index] = true;
                    current_free_variable = current_free_variable->next;
                }
                backtrack();
                return;
            }
        }
    }

    // We no longer need the list of infeasibility reducing variables
    delete_all_infeasibility_reducing_variables();
    
    /*
    There is a feasible continuation of the partial solution. We decide which variable to branch on using Balas' test.
    By convention, infeasibility reductions are negative if setting a variable to 1 reduces the infeasibility. 
    We already computed all of the infeasibility reductions earlier so we don't need to do it here
    */
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
    if (current_objective_value < best_objective_value && head_row == NULL){
        best_objective_value = current_objective_value;
        update_incumbent_solution();
    }

    reset_infeasibility_reductions();
    reset_minimum_infeasibilities();
}

void read_problem_data(char *filename){
    FILE *input_file;
    bool input_dimensions_valid = true;
    input_file = fopen(filename, "r");
    if (fscanf(input_file, "%d", &number_of_variables) != 1){
        input_dimensions_valid = false;
        printf("Error reading number of variables.\n");
    }
    if (fscanf(input_file, "%d", &number_of_constraints) != 1){
        input_dimensions_valid = false;
        printf("Error reading number of constraints.\n");
    }
    if (fscanf(input_file, "%d", &number_of_nonzeros) != 1){
        input_dimensions_valid = false;
        printf("Error reading number of nonzeros in constarints.\n");
    }
    if (fscanf(input_file, "%d", &number_of_objective_nonzeros) != 1){
        input_dimensions_valid = false;
        printf("Error reading number of nonzeros in objective.\n");
    }
    constraint_matrix = (double*) malloc(sizeof(double) * number_of_nonzeros);
    row_starts = (int*) malloc(sizeof(int) * number_of_constraints + 1);
    objective_coefficients = (double*) malloc(sizeof(double) * number_of_objective_nonzeros);
    objective_indices = (int*) malloc(sizeof(int) * number_of_objective_nonzeros);
    variable_indices = (int*) malloc(sizeof(int) * number_of_nonzeros);
    right_hand_sides = (double*) malloc(sizeof(double) * number_of_constraints);
    incumbent_solution = (bool*) calloc(number_of_variables, sizeof(bool));
    current_solution = (bool*) calloc(number_of_variables, sizeof(bool));
    left_hand_sides = (double*) calloc(number_of_constraints, sizeof(double));
    rows_with_negative_coefficients = malloc(sizeof(struct row_with_negative_coefficients**) * number_of_variables);
    rows_with_nonzero_coefficients = malloc(sizeof(struct row_with_nonzero_coefficients**) * number_of_variables);

    int current_index; // used to convert from 1-indexing to 0-indexing
    for (int i = 0; i < number_of_objective_nonzeros; i++){
        if (fscanf(input_file, "%d %lf", &current_index, objective_coefficients + i) != 2){
            input_dimensions_valid = false;
            printf("Error reading objective coefficient %d.\n", i + 1);
        }
        objective_indices[i] = current_index - 1;
    }
    for (int i = 0; i < number_of_constraints; i++){
        if (fscanf(input_file, "%lf", right_hand_sides + i) != 1){
            input_dimensions_valid = false;
            printf("Error reading right-hand side %d.\n", i + 1);
        }
    }
    int current_row = 0;
    for (int i = 0; i < number_of_nonzeros; i++){
        int previous_row = current_row;
        if (fscanf(input_file, "%d %d %lf", &current_row, &current_index, constraint_matrix + i) != 3){
            input_dimensions_valid = false;
            printf("Error reading nonzero %d in constraint matrix.\n", i + 1);
        }
        variable_indices[i] = current_index - 1;
        if (current_row != previous_row){
            row_starts[previous_row] = i;
        }
    }
    row_starts[number_of_constraints] = number_of_nonzeros + 1;

    for (int variable_index = 0; variable_index < number_of_variables; variable_index++){
        for (int row_index = 0; row_index < number_of_constraints; row_index++){
            int row_start = row_starts[row_index];
            int next_row_start = row_starts[row_index + 1];
            for (int j = row_start; j < next_row_start; j++){
                if (variable_indices[j] == variable_index){
                    double coefficient = constraint_matrix[j];
                    insert_row_with_nonzero_coefficient_end(variable_index, row_index, coefficient);
                    if (coefficient < 0){
                        insert_row_with_negative_coefficient_end(variable_index, row_index, coefficient);
                    }
                }
            }
        }
    }

    fclose(input_file);
    if (!input_dimensions_valid){
        printf("Exiting.\n");
        exit(-1);
    }
}

void free_all_memory(){
    /*
    Frees all memory at the end of the algorithm.
    */
    delete_all_infeasibility_reducing_variables();
    delete_all_free_variables();
    delete_all_infeasible_rows();
    free(constraint_matrix);
    free(row_starts);
    free(objective_coefficients);
    free(objective_indices);
    free(variable_indices);
    free(right_hand_sides);
    free(incumbent_solution);
    free(current_solution);
    free(left_hand_sides);
    delete_all_row_lists();
    free(rows_with_nonzero_coefficients);
    free(rows_with_negative_coefficients);
}