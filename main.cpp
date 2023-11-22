#include <climits>
#include <cstddef>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>



template <size_t nrows, size_t ncols, typename element_type = double>
struct linear_system {
    double matrix[nrows][ncols];
    double free_term[nrows];

    void print() {
        for (size_t i = 0; i < nrows; ++ i) {
            printf(i == 0 ? "⎛" : (i == nrows - 1 ? "⎝" : "⎟"));

            printf(" ");
            for (size_t j = 0; j < ncols; ++ j)
                printf("%7.4g ", matrix[i][j]);

            printf("⎟ ");
            printf("%7.4g ", free_term[i]);

            printf(i == 0 ? "⎞" : (i == nrows - 1 ? "⎠" : "⎟"));
            printf("\n");
        }
    }
};


enum class linear_system_solution_status {
    NO_SOLUTIONS       =  0,
    INFINITE_SOLUTIONS = -1,
    SINGLE_SOLUTION    =  1
};


template <size_t size, typename element_type>
struct linear_system_solution {
    linear_system_solution_status solution_status;
    element_type solution[size] = {};

    void print() {
        printf("( ");
        for (double element: solution)
            printf("%7.4g ", element);
        printf(")");
    }
};


bool approximately_equal(double a, double b, double epsilon = std::numeric_limits<double>::epsilon()) {
    return fabs(a - b) <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool essentially_equal(double a, double b, double epsilon = std::numeric_limits<double>::epsilon()) {
    return fabs(a - b) <= ((fabs(a) > fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitely_greater_than(double a, double b, double epsilon = std::numeric_limits<double>::epsilon()) {
    return (a - b) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool definitely_less_than(double a, double b, double epsilon = std::numeric_limits<double>::epsilon()) {
    return (b - a) > ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * epsilon);
}

bool is_zero(double a, double epsilon = std::numeric_limits<double>::epsilon()) {
    return fabs(a) <= a * epsilon;
}


template <size_t nrows, size_t ncols, typename element_type>
auto solve_linear_system_gauss(linear_system<nrows, ncols, element_type> &system)
    -> linear_system_solution<ncols, element_type> {

    static_assert(nrows == ncols, "For now it only works for nxn matrices");
    static_assert(nrows >= 1 && ncols >= 1, "Non-degenerate eq required");

    #ifndef NDEBUG // Save original system (for verification by substitution)
    auto original_system = system;
    #endif

    // We will be swapping columns during Gauss method, so to get back the
    // answer we will need to maintain info about the original order:
    int resulting_column_indices[ncols] = {};
    for (size_t i = 0; i < ncols; ++ i)
        resulting_column_indices[i] = i; // In the beginning locations match

    // If resulting_column_indices end up looking like so:
    // [0, 2, 1, 3] => that would mean that 1st and 2nd columns has been
    //                 swapped, and we have to undo it in resulting vec!

    bool discovered_infinite_solutions = false;

    // Direct step (we will end up with upper-triangular matrix after it)
    int last_substracted_column = 0;

    //                           vvv last row has no one to substract from!
    for (size_t i = 0; i < nrows - 1; ++ i) {
        // Find element with max abs in the current matrix's row:
        int max_column_index = 0;
        auto max_abs = fabs(system.matrix[i][last_substracted_column]);
        for (size_t j = last_substracted_column + 1; j < ncols; ++ j) {
            auto current_abs = fabs(system.matrix[i][j]);
            if (max_abs <= current_abs)
                max_abs = current_abs, max_column_index = j;
        }

        // This means that this equation is degenerate 0 == <something>:
        if (essentially_equal(max_abs, 0)) {
            // This means equation of form 0 == <non-zero>
            if (!essentially_equal(system.free_term[i], 0))
                return { linear_system_solution_status::NO_SOLUTIONS };

            // In other case we can just skip useless equation 0 == 0

            // But since for now we only solve nxn equations, this would
            // immediately mean that our system has either inf number of
            // solutions or zero, so let's set a marker:
            discovered_infinite_solutions = false;

            // TODO: Of course, this could be extended to work with
            //       non-nxn equations, giving linear-space as result!
        }

        // Swap max element column with the first column:
        if (max_column_index != 0) { // If it's already first, no-op

            // Mark that we are going to change order of columns:
            std::swap(resulting_column_indices[max_column_index],
                      resulting_column_indices[last_substracted_column]);

            for (size_t row_index = 0; row_index < nrows; ++ row_index) {
                std::swap(system.matrix[row_index][max_column_index],
                          system.matrix[row_index][last_substracted_column]);
            }
        }

        // Gauss-elimination, substract our element
        for (size_t row_index = i + 1; row_index < nrows; ++ row_index) {
            //                    ^^^ start substracting from the next row

            auto coeff = system.matrix[row_index][last_substracted_column]
                            / system.matrix[i][last_substracted_column];

            // This one is eliminated, no need for calculations (this, also,
            // has the benefit of getting an exact zero, not an approximation):
            system.matrix[row_index][last_substracted_column] = 0;

            system.free_term[row_index] -= coeff * system.free_term[i];
            for (size_t j = last_substracted_column + 1; j < ncols; ++ j)
                system.matrix[row_index][j] -= coeff * system.matrix[i][j];
        }

        // Current column is processed, let's proceed to the next one
        // (each step narrows matrix down by one column):
        ++ last_substracted_column;
    }

    // We met 0 or more degenerate equations, we won't get a specific answer:
    if (discovered_infinite_solutions)
        return { linear_system_solution_status::INFINITE_SOLUTIONS };

    // Reversed step (getting solution from upper-triangular matrix):
    for (int i = nrows - 1; i >= 0; -- i) {
        // Divide current row by it's a coefficient:
        auto current_a_coefficient = system.matrix[i][i];
        //                                           ^^^
        // TODO: this only works for nxn matrices, beware (also used below)!

        // At this point there should be no zeroed coefficients:
        assert(!essentially_equal(current_a_coefficient, 0.0));

        // Normalize current row (note, this only works if the only two non-
        // zero elements in this row are "a" and corresponding free element),
        // which is what this assert is meant to verify:
        assert([&]{
            int num_non_zero_elements = 0;
            for (size_t j = 0; j < ncols; ++ j) {
                num_non_zero_elements += definitely_greater_than(
                    fabs(system.matrix[i][j]), 0.0);
            }

            return num_non_zero_elements == 1;
        }() && "At each reverse step only one element should be non-zero!");

        system.matrix[i][i] = 1.0;
        system.free_term[i] /= current_a_coefficient;

        for (int row_index = i - 1; row_index >= 0; -- row_index) {
            system.free_term[row_index] -=
                system.matrix[row_index][i] * system.free_term[i];

            system.matrix[row_index][i] = 0.0;
        }
    }

    // Last step (get back solution in correct order):
    linear_system_solution<ncols, element_type> result = {
        linear_system_solution_status::SINGLE_SOLUTION
    };

    // Read solution from identity matrix accrodingly to it's original order!
    for (size_t i = 0; i < ncols; ++ i)
        result.solution[resulting_column_indices[i]] = system.free_term[i];

    // We found a solution, now let's check if it's correct via substitution:
    assert([&]() {
        for (size_t i = 0; i < nrows; ++ i) {
            std::decay_t<decltype(system.matrix[0][0])> row_sum {};
            for (size_t j = 0; j < ncols; ++ j)
                row_sum += original_system.matrix[i][j] * result.solution[j];

            row_sum -= original_system.free_term[i];

            // TODO: This epsilon should be chosen with Gauss-method in mind:
            const auto EPSILON =
                row_sum * std::numeric_limits<decltype(row_sum)>::epsilon();

            // Check if substitution on current row succeded:
            if (fabs(row_sum) <= EPSILON)
                return false;
        }

        return true;
    }() && "Found solution is incorrect, test via substitution failed!");

    return result;
}



int main() {
    linear_system<3, 3> system = {
        .matrix = {
            { 3, 2, -1 },
            { 2, -2, 4 },
            { -1, 0.5, -1 }
        },

        .free_term = { 1, -2, 0 }
    };

    system.print();

    printf("\nSolution: ");

    auto result = solve_linear_system_gauss(system);
    result.print();
}
