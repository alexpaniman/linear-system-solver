#include <climits>
#include <cstddef>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <limits>
#include <concepts>
#include <sstream>
#include <string>
#include <sys/wait.h>
#include <cppgfplots.h>



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


template <typename type>
struct latexify {
    static_assert(false, "There's no latexify specialization for your type!");
};


template <>
struct latexify<double> {
    std::string operator()(double value) {

        std::stringstream ss;
        ss << std::setprecision(2) << value;

        std::string result = ss.str();

        auto e_position = result.find('e');
        if (e_position != std::string::npos) {
            if (result[e_position + 1] == '+')
                result.erase(result.begin() + (e_position + 1));

            result.erase(result.begin() + e_position);

            result.insert(e_position + 1, "\\cdot 10^{");
            result.append("}");
        }

        return result;
    }
};

decltype(auto) to_latex(auto &&value) {
    return latexify<std::decay_t<decltype(value)>>{}(std::forward<decltype(value)>(value));
}



enum class linear_system_op { NONE = 0, MULTIPLIED_SUBSTRACT, COLUMN_SWAP, ROW_DIVIDE };

enum class linear_system_op_kind { COLUMN, ROW };

template <size_t nrows, size_t ncols, typename element_type>
class linear_system_viz {
public:
    linear_system_viz(linear_system<nrows, ncols> system):
        underlying_system_(system), current_operation_(linear_system_op::NONE) {};

    std::string latex() {
        std::stringstream ss;
        ss << "\\adjustbox{stack=cc}{\\vbox to -1.5ex {} \\(\n";

        ss << "\\begin{pNiceArray}[create-medium-nodes, first-row, last-row, ";
        if (current_operation_ == linear_system_op::MULTIPLIED_SUBSTRACT ||
            current_operation_ == linear_system_op::ROW_DIVIDE)
            ss << "last-col, ";

        ss << "extra-left-margin=.5em, extra-right-margin=.5em]{";

        for (size_t i = 0; i < ncols; ++ i)
            ss << "r";

        ss << "@{\\hskip 1em}|r}\n";

        {
            for (size_t i = 0; i < nrows; ++ i)
                ss << column_names_[i] << " &";

            ss << "\\\\\n";
        }

        for (size_t i = 0; i < nrows; ++ i) {
            ss << to_latex(underlying_system_.matrix[i][0]) << " ";
            for (size_t j = 1; j < ncols; ++ j)
                ss << "& " << to_latex(underlying_system_.matrix[i][j]) << " ";

            ss << "& " << to_latex(underlying_system_.free_term[i]) << " ";

            if (current_operation_ == linear_system_op::MULTIPLIED_SUBSTRACT
                || current_operation_ == linear_system_op::ROW_DIVIDE) {

                if (main_index != (int) i + 1)
                    ss << "& ";
                else {
                    double coefficient = coefficient_;
                    if (current_operation_ == linear_system_op::MULTIPLIED_SUBSTRACT) {
                        coefficient *= -1;
                        ss << "& \\times";
                    }

                    if (current_operation_ == linear_system_op::ROW_DIVIDE)
                        ss << " & \\div";

                    if (coefficient < 0)
                        ss << "(";

                    ss << to_latex(coefficient);

                    if (coefficient < 0)
                        ss << ")";
                }
            }


            ss << "\\\\\n";
        }

        ss << "\\\\\n";

        if (current_operation_ == linear_system_op::MULTIPLIED_SUBSTRACT ||
            current_operation_ == linear_system_op::ROW_DIVIDE) {

            glm::ivec2 sel0[] = {
                { main_index,           1 },
                { main_index,   ncols + 1 }
            };

            glm::ivec2 sel1[] = {
                { target_index,         1 },
                { target_index, ncols + 1 }
            };

            ss << "\\CodeAfter\n";

            ss << "\\tikz \\node [rectangle, draw=black, fit=";
            ss << "(" << sel0[0].x << "-" << sel0[0].y << "-medium)";
            ss << " ";
            ss << "(" << sel0[1].x << "-" << sel0[1].y << "-medium)";
            ss << "] {};\n";

            if (current_operation_ == linear_system_op::MULTIPLIED_SUBSTRACT) {
                ss << "\\tikz \\node [rectangle, draw=black, fit=";
                ss << "(" << sel1[0].x << "-" << sel1[0].y << "-medium)";
                ss << " ";
                ss << "(" << sel1[1].x << "-" << sel1[1].y << "-medium)";
                ss << "] {};\n";

                ss << "\\tikz \\node [at=";
                ss << "(" << main_index << ".5-|" << (ncols + 2) << ")";
                ss << ", right=1em] (main-row-proxy) {};\n";

                if (main_index < target_index)
                    ss << "\\tikz \\node [at=(main-row-proxy), below=.1ex] (main-row) {};";
                else {
                    ss << "\\tikz \\node [at=(main-row-proxy)] (main-row) {};";
                }

                ss << "\\tikz \\node [at=";
                ss << "(" << target_index << ".5-|" << (ncols + 2) << ")";
                ss << ", right=.6em] (target-row-proxy) {};\n";

                if (main_index < target_index)
                    ss << "\\tikz \\node [at=(target-row-proxy),  below=.15ex] (target-row) {};";
                else {
                    ss << "\\tikz \\node [at=(target-row-proxy), above=.1ex] (target-row) {};";
                }

                ss << "\\tikz \\path [thick, ->] ";
//              ss << "(" << main_index <<      "-" << (ncols + 2) << "-medium)";
                ss << "(main-row) edge[pos=0.5, \"\\texttt{+}\", "
                   << (main_index < target_index ? "bend left" : "bend right")
                   << "] (target-row);\n";
            }
        }

        if (current_operation_ == linear_system_op::COLUMN_SWAP) {

            // TODO: cleanup, node generation already happends in other
            //       two operations no need to create a completely new
            //       case for them!
            ss << "\\CodeAfter\n";

            int first_index  = std::min(main_index, target_index);
            int second_index = std::max(main_index, target_index);

            ss << "\\tikz \\node [rectangle, draw=black, fit=";
            ss << "(" <<     1 << "-" << first_index << "-medium)";
            ss << " ";
            ss << "(" << ncols << "-" << first_index << "-medium)";
            ss << "] {};\n";

            ss << "\\tikz \\node [rectangle, draw=black, fit=";
            ss << "(" <<     1 << "-" << second_index << "-medium)";
            ss << " ";
            ss << "(" << ncols << "-" << second_index << "-medium)";
            ss << "] {};\n";

            ss << "\\tikz \\path [thick, <->] ";
            ss << "(0-" << main_index << "-medium)";
            ss << " edge[bend right] (0-" << target_index << "-medium);\n";
        }

        ss << "\\end{pNiceArray} \\mkern1mu";
        ss << "\\) \\\\ \\vbox to -1.5ex {}}";


        return ss.str();
    }

    linear_system<nrows, ncols, element_type> underlying_system_;

    linear_system_op current_operation_ = linear_system_op::NONE;

    double coefficient_ = 0.0;
    linear_system_op_kind operation_kind_ = linear_system_op_kind::ROW;

    std::string column_names_[ncols];

    int main_index = 0, target_index = 0;
};



template <size_t nrows, size_t ncols, typename element_type>
class linear_system_transform_path {
public:
    std::vector<linear_system_viz<nrows, ncols, element_type>> transforms_;

    void generate_report(std::string output_file_name) {
        std::filesystem::path path = output_file_name;
        if (!path.parent_path().empty() && !std::filesystem::exists(path.parent_path()))
            execute("mkdir", ".", path.parent_path().c_str());

        char temp_directory_name[] = "/tmp/linear-system-visualization-XXXXXX";
        mkdtemp(temp_directory_name);

        std::string image_tex = temp_directory_name;
        image_tex += "/standalone-image.tex";

        FILE *image_tex_stream = fopen(image_tex.c_str(), "w");
        fprintf(image_tex_stream, 
            "\\documentclass[12pt]{extarticle}"                                  "\n"
                                                                                        "\n"
            "\\tolerance=10000"                                                         "\n"
                                                                                        "\n"
            "\\usepackage[paperwidth=60in,paperheight=100in,left=3cm,right=3cm,top=1cm,bottom=1cm]{geometry}"             "\n"
                                                                                        "\n"
            "\\usepackage{nicematrix}"                                                  "\n"
            "\\usepackage{adjustbox}"                                                   "\n"
                                                                                        "\n"
            "\\usepackage{tikz}"                                                        "\n"
            "\\usetikzlibrary{fit,shapes.geometric,quotes}"                             "\n"
R"(
\tikzset{C/.style={circle, draw}}
\tikzset{every edge quotes/.style = {circle, fill=white, inner sep=.04em, anchor=center}}
)"
                                                                                        "\n"
            "\\linespread{1.5}"                                                         "\n"
                                                                                        "\n"
            "\\begin{document}"                                                         "\n"
            "\\setlength\\arrayrulewidth{0.6pt}"                                        "\n"
        //  "\\setlength\\arraycolsep{.5em}"                                            "\n"
                                                                                        "\n"
        );

        fprintf(image_tex_stream, "\\noindent");
        for (size_t i = 0; i < transforms_.size(); ++i) {
            fprintf(image_tex_stream, "%s", transforms_[i].latex().c_str());
            if (i != transforms_.size() - 1)
                fprintf(image_tex_stream, "~\\(\\boldsymbol{\\sim}\\)\n");
        }

        fprintf(image_tex_stream,
                "\\end{document}"                                                       "\n"
        );

        fclose(image_tex_stream), image_tex_stream = NULL;

        execute("latexmk", temp_directory_name, "-pdf",
                "-pdflatex=lualatex -interaction=nonstopmode %O %S", image_tex.c_str());

        std::string image_pdf = temp_directory_name;
        image_pdf += "/standalone-image.pdf";

        execute("mv", ".", image_pdf.c_str(), output_file_name.c_str());
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
auto solve_linear_system_gauss_impl(linear_system<nrows, ncols, element_type> &system,
                                    auto &&report_matrix_transformation)
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

    for (size_t i = 0; i < nrows; ++ i) {
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
            {
                linear_system_viz<nrows, ncols, element_type> current_viz{system};
                current_viz.current_operation_ = linear_system_op::COLUMN_SWAP;
                current_viz.  main_index = max_column_index + 1;
                current_viz.target_index = last_substracted_column + 1;

                report_matrix_transformation(current_viz);
            }

            // Mark that we are going to change order of columns:
            std::swap(resulting_column_indices[max_column_index],
                      resulting_column_indices[last_substracted_column]);

            for (size_t row_index = 0; row_index < nrows; ++ row_index) {
                std::swap(system.matrix[row_index][max_column_index],
                          system.matrix[row_index][last_substracted_column]);
            }
        }

        {
            linear_system_viz<nrows, ncols, element_type> current_viz{system};
            current_viz.current_operation_ = linear_system_op::ROW_DIVIDE;
            current_viz.  main_index = last_substracted_column + 1;
            current_viz.coefficient_ = system.matrix[i][last_substracted_column];
            report_matrix_transformation(current_viz);
        }

        // Normalize current row:
        for (size_t j = last_substracted_column + 1; j < ncols; ++ j)
            system.matrix[i][j] /= system.matrix[i][last_substracted_column];

        system.free_term[i] /= system.matrix[i][last_substracted_column];
        system.matrix[i][last_substracted_column] = 1.0;

        // Gauss-elimination, substract our element
        for (size_t row_index = i + 1; row_index < nrows; ++ row_index) {
            //                    ^^^ start substracting from the next row

            auto coeff = system.matrix[row_index][last_substracted_column];

            {
                linear_system_viz<nrows, ncols, element_type> current_viz{system};
                current_viz.current_operation_ = linear_system_op::MULTIPLIED_SUBSTRACT;
                current_viz.  main_index =         i + 1;
                current_viz.target_index = row_index + 1;
                current_viz.coefficient_ = coeff;

                report_matrix_transformation(current_viz);
            }

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

        for (int row_index = i - 1; row_index >= 0; -- row_index) {
            {
                linear_system_viz<nrows, ncols, element_type> current_viz{system};
                current_viz.current_operation_ = linear_system_op::MULTIPLIED_SUBSTRACT;
                current_viz.  main_index = i + 1;
                current_viz.target_index = row_index + 1;
                current_viz.coefficient_ = system.matrix[row_index][i];

                report_matrix_transformation(current_viz);
            }

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

            // TODO: This epsilon should be chosen with Gauss-method in mind
            //       (now It's just a trial-and-error-selected arbitrary value):
            // Check if substitution on current row succeded:
            if (fabs(row_sum) >= 1e-1)
                return false;
        }

        return true;
    }() && "Found solution is incorrect, test via substitution failed!");


    {
        linear_system_viz<nrows, ncols, element_type> current_viz{system};
        for (size_t i = 0; i < ncols; ++ i)
            current_viz.column_names_[i] = "x^{" + std::to_string(resulting_column_indices[i] + 1) + "}";

        report_matrix_transformation(current_viz);
    }

    return result;
}

template <size_t nrows, size_t ncols, typename element_type>
decltype(auto) solve_linear_system_gauss(linear_system<nrows, ncols, element_type> &system) {
    return solve_linear_system_gauss_impl(system,
            []([[maybe_unused]] linear_system_viz<nrows, ncols, element_type> viz){}
    );
}

template <size_t nrows, size_t ncols, typename element_type>
auto solve_linear_system_gauss(linear_system<nrows, ncols, element_type> &system,
                                         std::string output_file_name) {

    linear_system_transform_path<nrows, ncols, element_type> path {};
    auto result = solve_linear_system_gauss_impl(system,
            [&](linear_system_viz<nrows, ncols, element_type> viz){
                path.transforms_.emplace_back(std::move(viz));
            }
    );

    path.generate_report(output_file_name);

    return result;
}



int main() {
    linear_system<10, 10> system = {};
    for (size_t i = 0; i < 10; ++ i) {
        for (size_t j = 0; j < 10; ++ j)
            system.matrix[i][j] = rand() + rand() / (double)INT_MAX;

        system.free_term[i] = rand() + rand() / (double)INT_MAX;
    }

    std::cout << "\nSolution: ";
    std::cout.flush();

    auto result = solve_linear_system_gauss(system, "gauss-report.pdf");
    result.print();
}
