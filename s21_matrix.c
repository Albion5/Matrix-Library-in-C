#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
  int status = OK;
  matrix_t error = {NULL, 0, 0};
  if (rows > 0 && columns > 0 && result != NULL) {
    result->rows = rows;
    result->columns = columns;
    result->matrix = calloc(rows, sizeof(double *));
    for (int i = 0; i < rows; i++) {
      result->matrix[i] = calloc(columns, sizeof(double));
    }
  } else {
    if (result) *result = error;
    status = INCORRECT_M;
  }
  return status;
}

int is_valid_matrix(matrix_t *result) {
  return result != NULL && result->matrix != NULL && result->rows > 0 &&
         result->columns > 0;
}

void s21_remove_matrix(matrix_t *A) {
  if (A != NULL && A->matrix != NULL) {
    for (int i = 0; i < A->rows; i++) {
      if (A->matrix[i] != NULL) {
        free(A->matrix[i]);
      }
    }
    free(A->matrix);
  }
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
  int status = is_valid_matrix(A) && is_valid_matrix(B) && is_eq_size(A, B);
  if (status) {
    for (int row = 0; row < A->rows; row++) {
      for (int column = 0; column < A->columns; column++) {
        if (fabs(A->matrix[row][column] - B->matrix[row][column]) >= 1e-7) {
          status = FAILURE;
        }
      }
    }
  }
  return status;
}

int is_eq_size(matrix_t *A, matrix_t *B) {
  return (A->columns == B->columns && A->rows == B->rows);
}

int s21_sum_common(matrix_t *A, matrix_t *B, matrix_t *result, int k) {
  int status = OK;
  matrix_t error = {NULL, 0, 0};
  if (is_valid_matrix(A) && is_valid_matrix(B) && result != NULL) {
    if (is_eq_size(A, B)) {
      s21_create_matrix(A->rows, A->columns, result);
      for (int row = 0; row < A->rows; row++) {
        for (int column = 0; column < A->columns; column++) {
          result->matrix[row][column] =
              A->matrix[row][column] + B->matrix[row][column] * k;
        }
      }
    } else {
      status = CALC_ERROR;
      *result = error;
    }
  } else {
    status = INCORRECT_M;
    if (result) *result = error;
  }
  return status;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = s21_sum_common(A, B, result, 1);
  return status;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = s21_sum_common(A, B, result, -1);
  return status;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
  int status = OK;
  matrix_t error = {NULL, 0, 0};
  if (is_valid_matrix(A) && result != NULL) {
    s21_create_matrix(A->rows, A->columns, result);
    for (int row = 0; row < A->rows; row++) {
      for (int column = 0; column < A->columns; column++) {
        result->matrix[row][column] = A->matrix[row][column] * number;
      }
    }
  } else {
    status = INCORRECT_M;
    if (result) *result = error;
  }
  return status;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
  int status = OK;
  matrix_t error = {NULL, 0, 0};
  if (is_valid_matrix(A) && is_valid_matrix(B) && result != NULL) {
    if (A->columns == B->rows) {
      s21_create_matrix(A->rows, B->columns, result);
      for (int row = 0; row < result->rows; row++) {
        for (int column = 0; column < result->columns; column++) {
          result->matrix[row][column] =
              multiply_row_by_column(A, B, A->columns, row, column);
        }
      }
    } else {
      status = CALC_ERROR;
      *result = error;
    }
  } else {
    status = INCORRECT_M;
    if (result) *result = error;
  }
  return status;
}

double multiply_row_by_column(matrix_t *A, matrix_t *B, int times, int row,
                              int column) {
  double n = 0;
  for (int i = 0; i < times; i++) {
    n += A->matrix[row][i] * B->matrix[i][column];
  }
  return n;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
  int status = OK;
  matrix_t error = {NULL, 0, 0};
  if (is_valid_matrix(A) && result != NULL) {
    s21_create_matrix(A->columns, A->rows, result);
    for (int row = 0; row < result->rows; row++) {
      for (int column = 0; column < result->columns; column++) {
        result->matrix[row][column] = A->matrix[column][row];
      }
    }
  } else {
    status = INCORRECT_M;
    if (result) *result = error;
  }
  return status;
}

int get_sign(int row) { return (row % 2) ? -1 : 1; }

int s21_calc_complements(matrix_t *A, matrix_t *result) {
  int status = OK;
  matrix_t error = {NULL, 0, 0};
  if (is_valid_matrix(A) && result != NULL) {
    if (A->columns == A->rows) {
      matrix_t minor_matrix;
      int size = A->rows;
      s21_create_matrix(size, size, result);
      s21_create_matrix(size - 1, size - 1, &minor_matrix);
      for (int row = 0; row < result->rows; row++) {
        for (int column = 0, sign = get_sign(row); column < result->columns;
             column++, sign = -sign) {
          get_minor_matrix(A, &minor_matrix, row, column);
          result->matrix[row][column] = sign * get_determinant(&minor_matrix);
        }
      }
      s21_remove_matrix(&minor_matrix);
    } else {
      status = CALC_ERROR;
      *result = error;
    }
  } else {
    status = INCORRECT_M;
    if (result) *result = error;
  }
  return status;
}

void get_minor_matrix(matrix_t *A, matrix_t *minor_matrix, int i, int j) {
  for (int row = 0, m_row = 0; row < A->rows; row++) {
    for (int column = 0, m_column = 0; column < A->columns; column++) {
      if (row != i && column != j) {
        minor_matrix->matrix[m_row][m_column++] = A->matrix[row][column];
      }
    }
    if (row != i) m_row++;
  }
}

double get_determinant(matrix_t *A) {
  double det = 0;
  int size = A->rows;
  if (size == 1) {
    det = A->matrix[0][0];
  } else {
    matrix_t minor_matrix;
    s21_create_matrix(size - 1, size - 1, &minor_matrix);
    double *first_row = A->matrix[0];
    for (int sign = 1, column = 0; column < size; column++, sign = -sign) {
      get_minor_matrix(A, &minor_matrix, 0, column);
      det += sign * first_row[column] * get_determinant(&minor_matrix);
    }
    s21_remove_matrix(&minor_matrix);
  }
  return det;
}

int s21_determinant(matrix_t *A, double *result) {
  int status = OK;
  if (is_valid_matrix(A) && result != NULL) {
    if (A->columns == A->rows) {
      *result = get_determinant(A);
    } else {
      status = CALC_ERROR;
    }
  } else {
    status = INCORRECT_M;
  }
  return status;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
  int status = OK;
  matrix_t error = {NULL, 0, 0};
  if (is_valid_matrix(A) && result != NULL) {
    if (A->columns == A->rows) {
      double det = get_determinant(A);
      if (det != 0) {
        matrix_t cofactor_matrix, transpose;
        s21_calc_complements(A, &cofactor_matrix);
        s21_transpose(&cofactor_matrix, &transpose);
        s21_mult_number(&transpose, 1 / det, result);
        s21_remove_matrix(&cofactor_matrix);
        s21_remove_matrix(&transpose);
      } else {
        status = CALC_ERROR;
        *result = error;
      }
    } else {
      status = CALC_ERROR;
      *result = error;
    }
  } else {
    status = INCORRECT_M;
    if (result) *result = error;
  }
  return status;
}

void create_test_matrix(const char *filename, int i, int j, matrix_t *A) {
  char path[100] = "tests/data/";
  strcat(path, filename);
  FILE *f = fopen(path, "r");
  s21_create_matrix(i, j, A);
  for (int row = 0; row < i; row++) {
    for (int column = 0; column < j; column++) {
      fscanf(f, "%lf", A->matrix[row] + column);
    }
  }
  fclose(f);
}
