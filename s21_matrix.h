#ifndef C6_S21_MATRIX_S21_MATRIX_H_
#define C6_S21_MATRIX_S21_MATRIX_H_
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SUCCESS 1
#define FAILURE 0
typedef struct matrix_struct {
  double **matrix;
  int rows;
  int columns;
} matrix_t;

enum { OK, INCORRECT_M, CALC_ERROR };
// All operations (except matrix comparison) should return the resulting code:
// 0 - OK
// 1 - Error, incorrect matrix
// 2 - Calculation error (mismatched matrix sizes; matrix for which calculations
// cannot be performed, etc.)

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);
int s21_eq_matrix(matrix_t *A, matrix_t *B);
int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_transpose(matrix_t *A, matrix_t *result);
int s21_inverse_matrix(matrix_t *A, matrix_t *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);
int s21_determinant(matrix_t *A, double *result);

// Helper functions
int is_valid_matrix(matrix_t *result);
int is_eq_size(matrix_t *A, matrix_t *B);
int s21_sum_common(matrix_t *A, matrix_t *B, matrix_t *result, int k);
double multiply_row_by_column(matrix_t *A, matrix_t *B, int times, int row,
                              int column);
double get_determinant(matrix_t *A);
void get_minor_matrix(matrix_t *A, matrix_t *minor_matrix, int i, int j);
void create_test_matrix(const char *filename, int i, int j, matrix_t *A);
#endif  // C6_S21_MATRIX_S21_MATRIX_H_