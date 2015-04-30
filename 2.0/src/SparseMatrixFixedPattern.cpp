#include "SparseMatrixFixedPattern.h"

SparseMatrixFixedPattern::SparseMatrixFixedPattern() {
    num_rows = 0;
    num_nonzero_elements = 0;
    entry = NULL;
    key = NULL;
    source_list = NULL;
    diagonal = NULL;
}

SparseMatrixFixedPattern::~SparseMatrixFixedPattern() {
    delete[] entry;
    delete[] key;
    delete[] source_list;
    delete[] diagonal;
    entry = NULL;
    key = NULL;
    source_list = NULL;
    diagonal = NULL;
    num_rows = 0;
    num_nonzero_elements = 0;
}

void SparseMatrixFixedPattern::init(int num_rows, int num_nonzero_elements, sparse_entry *entry, int *key, sparse_entry_sources *source_list) {
    this->num_rows = num_rows;
    this->num_nonzero_elements = num_nonzero_elements;
    this->entry = entry;
    this->key = key;
    this->source_list = source_list;

    // Work out which elements are on the diagonal
    diagonal = new scalar*[num_rows];
    for (int i = 0; i < num_rows; i++) {
        for (int j = key[i]; j < key[i + 1]; j++) {
            if (i == entry[j].column_index) {
                diagonal[i] = &(entry[j].val);
            }
        }
    }
}

/* Reconstruct the matrix by adding up all the contributions from the sources stored in the source list */
void SparseMatrixFixedPattern::build() {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i)
#endif
    for (i = 0; i < num_nonzero_elements; i++) {
        entry[i].val = source_list[i].sum_all_sources();
    }
}

/* Applies this matrix to the given vector 'in', writing the result to 'result' */
void SparseMatrixFixedPattern::apply(scalar *in, scalar *result) {
    for (int i = 0; i < num_rows; i++) {
        result[i] = 0;
        for (int j = key[i]; j < key[i + 1]; j++) {
            result[i] += entry[j].val * in[entry[j].column_index];
        }
    }
}

/* Applies this matrix to the given vector 'in', writing the result to 'result'. 'in' is made of 'vector3's */
void SparseMatrixFixedPattern::apply(vector3 *in, vector3 *result) {
    int i, j;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i, j) shared(result, in)
#endif
    for (i = 0; i < num_rows / 3; i++) {
        result[i].x = 0;
        result[i].y = 0;
        result[i].z = 0;
        for (j = key[3 * i]; j < key[(3 * i) + 1]; j++) {
            if (entry[j].column_index % 3 == 0) {
                result[i].x += entry[j].val * in[entry[j].column_index / 3].x;
            } else if (entry[j].column_index % 3 == 1) {
                result[i].x += entry[j].val * in[entry[j].column_index / 3].y;
            } else if (entry[j].column_index % 3 == 2) {
                result[i].x += entry[j].val * in[entry[j].column_index / 3].z;
            }
        }
        for (j = key[(3 * i) + 1]; j < key[(3 * i) + 2]; j++) {
            if (entry[j].column_index % 3 == 0) {
                result[i].y += entry[j].val * in[entry[j].column_index / 3].x;
            } else if (entry[j].column_index % 3 == 1) {
                result[i].y += entry[j].val * in[entry[j].column_index / 3].y;
            } else if (entry[j].column_index % 3 == 2) {
                result[i].y += entry[j].val * in[entry[j].column_index / 3].z;
            }
        }
        for (j = key[(3 * i) + 2]; j < key[(3 * i) + 3]; j++) {
            if (entry[j].column_index % 3 == 0) {
                result[i].z += entry[j].val * in[entry[j].column_index / 3].x;
            } else if (entry[j].column_index % 3 == 1) {
                result[i].z += entry[j].val * in[entry[j].column_index / 3].y;
            } else if (entry[j].column_index % 3 == 2) {
                result[i].z += entry[j].val * in[entry[j].column_index / 3].z;
            }
        }
    }
}

void SparseMatrixFixedPattern::calc_inverse_diagonal(scalar *inv_D) {
    int i;
#ifdef FFEA_PARALLEL_WITHIN_BLOB
#pragma omp parallel for default(none) private(i) shared(inv_D)
#endif
    for (i = 0; i < num_rows; i++) {
        inv_D[i] = 1.0 / (*(diagonal[i]));
    }
}

void SparseMatrixFixedPattern::print() {
    for (int i = 0; i < num_nonzero_elements; i++) {
        printf("[%d %e]\n", entry[i].column_index, entry[i].val);
    }
}

void SparseMatrixFixedPattern::print_dense() {
    int i, j, l;
    for (i = 0; i < num_rows; ++i) {
        l = 0;
        for (j = 0; j < num_rows; ++j) {
            if (j == entry[key[i] + l].column_index && key[i] + l < key[i + 1]) {
                printf("%e,", entry[key[i] + l].val);
                l++;
            } else {
                printf("%e,", 0.0);
            }
        }
        printf("\n");
    }
}

/* Prints dense matrix out to file for analysis. I suggest only letting this function run once (step = 1?) */
void SparseMatrixFixedPattern::print_dense_to_file(vector3 *a) {
    FILE *fout, *fout2;
    fout = fopen("dense_matrix.csv", "w");
    fout2 = fopen("force.csv", "w");
    int i, j, l;
    for (i = 0; i < num_rows; ++i) {
        l = 0;
        for (j = 0; j < num_rows; ++j) {
            if (j == entry[key[i] + l].column_index && key[i] + l < key[i + 1]) {
                fprintf(fout, "%e,", entry[key[i] + l].val);
                l++;
            } else {
                fprintf(fout, "%e,", 0.0);
            }
        }
        fprintf(fout, "\n");
    }
    for (i = 0; i < num_rows / 3; ++i) {
        fprintf(fout2, "%e\n%e\n%e\n", a[i].x, a[i].y, a[i].z);
    }
    fclose(fout);
    fclose(fout2);
}

void SparseMatrixFixedPattern::print_row_column() {
    FILE *fout;
    fout = fopen("row_column.csv", "w");
    for (int i = 0; i < num_rows; ++i) {
        for (int j = key[i]; j < key[i + 1]; ++j) {
            fprintf(fout, "%d,%d\n", entry[j].column_index, i);
        }
    }
    fclose(fout);
}

void SparseMatrixFixedPattern::check_symmetry() {
    int row;

    for (int i = 0; i < num_rows / 3; i++) {
        for (int j = key[i]; j < key[i + 1]; j++) {
            printf("%e ", entry[j].val);
            row = entry[j].column_index;
            for (int k = key[row]; k < key[row + 1]; ++k) {
                if (entry[k].column_index == i) {
                    printf("%e %e\n", entry[k].val, entry[k].val - entry[j].val);
                }
            }
        }
    }
}

void SparseMatrixFixedPattern::am_i_diagonally_dominant() {
    int i, j;
    scalar sum;
    FILE *fout;
    fout = fopen("diag_domin.csv", "w");
    for (i = 0; i < num_rows; ++i) {
        sum = 0.0;
        for (j = key[i]; j < key[i + 1]; ++j) {
            if (entry[j].column_index == i) {
                continue;
            }
            sum += fabs(entry[j].val);
        }
        fprintf(fout, "%d,%f\n", i, sum / fabs(*diagonal[i]));
    }
    fclose(fout);
}

// Works in progress

void SparseMatrixFixedPattern::cholesky_decompose(scalar* L) {
    int i, j, k = 0;

    // Recreating dense matrix
    for (i = 0; i < num_rows; ++i) {
        for (j = 0; j < num_rows; ++j) {
            if (j == entry[k].column_index) {
                L[i * num_rows + j] = entry[k].val;
                k++;
            } else {
                L[i * num_rows + j] = 0.0;
            }
        }
    }

    // Decomposition Loop
    for (i = 0; i < num_rows - 1; ++i) {
        L[i * num_rows + i] = sqrt(L[i * num_rows + i]);
        for (j = i + 1; j < num_rows; ++j) {
            L[j * num_rows + i] *= 1.0 / L[i * num_rows + i];
        }
        for (j = i + 1; j < num_rows; ++j) {
            for (k = i + 1; k < num_rows; ++k) {
                L[j * num_rows + k] -= L[j * num_rows + i] * L[k * num_rows + i];
            }
        }
    }

    // Zeroing upper triangle
    for (i = 1; i < num_rows; ++i) {
        for (j = 0; j < i; ++j) {
            L[j * num_rows + i] = 0.0;
        }
    }
}

void SparseMatrixFixedPattern::forwardbacksub(vector3* f) {
    int i, j, temp2;
    scalar temp;

    // Forward First
    for (i = 0; i < num_rows / 3; ++i) {
        temp = 0.0;
        for (j = key[3 * i]; j < key[3 * i + 1] - 1; ++j) {
            temp2 = entry[j].column_index % 3;
            if (temp2 == 0) {
                temp += entry[j].val * f[entry[j].column_index / 3].x;
            } else if (temp2 == 1) {
                temp += entry[j].val * f[entry[j].column_index / 3].y;
            } else {
                temp += entry[j].val * f[entry[j].column_index / 3].z;
            }
        }
        f[i].x = (f[i].x - temp) / *diagonal[3 * i];

        temp = 0.0;
        for (j = key[3 * i + 1]; j < key[3 * i + 2] - 1; ++j) {
            temp2 = entry[j].column_index % 3;
            if (temp2 == 0) {
                temp += entry[j].val * f[entry[j].column_index / 3].x;
            } else if (temp2 == 1) {
                temp += entry[j].val * f[entry[j].column_index / 3].y;
            } else {
                temp += entry[j].val * f[entry[j].column_index / 3].z;
            }
        }
        f[i].y = (f[i].y - temp) / *diagonal[3 * i + 1];

        temp = 0.0;
        for (j = key[3 * i + 2]; j < key[3 * i + 3] - 1; ++j) {
            temp2 = entry[j].column_index % 3;
            if (temp2 == 0) {
                temp += entry[j].val * f[entry[j].column_index / 3].x;
            } else if (temp2 == 1) {
                temp += entry[j].val * f[entry[j].column_index / 3].y;
            } else {
                temp += entry[j].val * f[entry[j].column_index / 3].z;
            }
        }
        f[i].z = (f[i].z - temp) / *diagonal[3 * i + 2];
    }

    // Now back
    for (i = num_rows / 3 - 1; i >= 0; --i) {
        temp = 0.0;
        for (j = key[3 * i + 3] - 1; j >= key[3 * i + 2]; --j) {
            temp2 = entry[j].column_index % 3;
            if (temp2 == 0) {
                temp += entry[j].val * f[entry[j].column_index / 3].x;
            } else if (temp2 == 1) {
                temp += entry[j].val * f[entry[j].column_index / 3].y;
            } else {
                temp += entry[j].val * f[entry[j].column_index / 3].z;
            }
        }
        f[i].z = (f[i].z - temp) / *diagonal[3 * i + 2];

        temp = 0.0;
        for (j = key[3 * i + 2] - 1; j >= key[3 * i + 1]; --j) {
            temp2 = entry[j].column_index % 3;
            if (temp2 == 0) {
                temp += entry[j].val * f[entry[j].column_index / 3].x;
            } else if (temp2 == 1) {
                temp += entry[j].val * f[entry[j].column_index / 3].y;
            } else {
                temp += entry[j].val * f[entry[j].column_index / 3].z;
            }
        }
        f[i].y = (f[i].y - temp) / *diagonal[3 * i + 1];

        temp = 0.0;
        for (j = key[3 * i + 1] - 1; j >= key[3 * i]; --j) {
            temp2 = entry[j].column_index % 3;
            if (temp2 == 0) {
                temp += entry[j].val * f[entry[j].column_index / 3].x;
            } else if (temp2 == 1) {
                temp += entry[j].val * f[entry[j].column_index / 3].y;
            } else {
                temp += entry[j].val * f[entry[j].column_index / 3].z;
            }
        }
        f[i].x = (f[i].x - temp) / *diagonal[3 * i];
    }
}

