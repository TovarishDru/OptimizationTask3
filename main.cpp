#include <string>
#include <cmath>
#include <vector>
#include <iostream>
#include<algorithm>
using namespace std;
using ll = long long;
using ull = unsigned long long;
using db = double;
using uint = unsigned int;


class Matrix {
protected:
    int n;
    int m;
    unique_ptr<unique_ptr<double[]>[]> matrix;
    friend ostream& operator<<(ostream& os, const Matrix& matrix) {
        string output = "";
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                string num = to_string(round(matrix.matrix[i][j] * matrix.roundValue) / matrix.roundValue);
                for (int k = 0; k < num.size(); k++) {
                    output += num[k];
                    if (num[k] == '.') {
                        for (int f = 0; f < log10(matrix.roundValue); f++) {
                            output += num[k + f + 1];
                        }
                        break;
                    }
                }
                if (j < matrix.m - 1) {
                    output += " ";
                }
            }
            output += "\n";
        }
        return os << output;
    }
    friend istream& operator>>(istream& is, const Matrix& matrix) {
        for (int i = 0; i < matrix.n; i++) {
            for (int j = 0; j < matrix.m; j++) {
                is >> matrix.matrix[i][j];
            }
        }
        return is;
    }
public:
    int roundValue = 1e6;
    Matrix operator+(const Matrix& matrix) const {
        if (this->n == matrix.n and this->m == matrix.m) {
            Matrix res(this->n, this->m);
            for (int i = 0; i < this->n; i++) {
                for (int j = 0; j < this->m; j++) {
                    res.matrix[i][j] = this->matrix[i][j] + matrix.matrix[i][j];
                }
            }
            return res;
        }
        else {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
    }
    Matrix* operator+(const Matrix* matrix) const {
        return new Matrix(*this + *matrix);
    }
    Matrix operator-(const Matrix& matrix) const {
        if (this->n == matrix.n and this->m == matrix.m) {
            Matrix res(this->n, this->m);
            for (int i = 0; i < this->n; i++) {
                for (int j = 0; j < this->m; j++) {
                    res.matrix[i][j] = this->matrix[i][j] - matrix.matrix[i][j];
                }
            }
            return res;
        }
        else {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
    }
    Matrix* operator-(Matrix* matrix) const {
        return new Matrix(*this - *matrix);
    }
    Matrix operator*(const Matrix& matrix) const {
        if (this->m == matrix.n) {
            Matrix res(this->n, matrix.m);
            for (int i = 0; i < this->n; i++) {
                for (int j = 0; j < matrix.m; j++) {
                    for (int k = 0; k < this->m; k++) {
                        res.matrix[i][j] += this->matrix[i][k] * matrix.matrix[k][j];
                    }
                }
            }
            return res;
        }
        else {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
    }
    Matrix operator*(const double& constant) {
        Matrix res(this->n, this->m);
        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < this->m; j++) {
                res.matrix[i][j] = this->matrix[i][j] * constant;
            }
        }
        return res;
    }
    Matrix* operator*(const Matrix* matrix) const {
        return new Matrix(*this * *matrix);
    }
    Matrix operator=(const Matrix& matrix) {
        this->n = matrix.n;
        this->m = matrix.m;
        this->matrix.release();
        this->matrix = unique_ptr<unique_ptr<double[]>[]>(new unique_ptr<double[]>[this->n]);
        for (int i = 0; i < matrix.n; i++) {
            this->matrix[i] = unique_ptr<double[]>(new double[this->m]);
            for (int j = 0; j < matrix.m; j++) {
                this->matrix[i][j] = matrix.matrix[i][j];
            }
        }
        return *this;
    }
    Matrix* operator=(const Matrix* matrix) {
        *this = *matrix;
        return this;
    }
    Matrix transpose() const {
        Matrix res(this->m, this->n);
        for (int i = 0; i < this->n; i++) {
            for (int j = 0; j < this->m; j++) {
                res.matrix[j][i] = this->matrix[i][j];
            }
        }
        return res;
    }
    int getN() const {
        return this->n;
    }
    int getM() const {
        return this->m;
    }
    double getElem(int i, int j) const {
        if (i >= this->n or j >= this->m) {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
        return this->matrix[i][j];
    }
    void setElem(int i, int j, double val) {
        this->matrix[i][j] = val;
    }
    Matrix() = delete;
    Matrix(int n, int m) {
        this->n = n;
        this->m = m;
        matrix = unique_ptr<unique_ptr<double[]>[]>(new unique_ptr<double[]>[n]);
        for (int i = 0; i < n; i++) {
            matrix[i] = unique_ptr<double[]>(new double[m]);
            for (int j = 0; j < m; j++) {
                matrix[i][j] = 0;
            }
        }
    }
    Matrix(const Matrix& matrix) {
        this->n = matrix.n;
        this->m = matrix.m;
        this->matrix = unique_ptr<unique_ptr<double[]>[]>(new unique_ptr<double[]>[this->n]);
        for (int i = 0; i < matrix.n; i++) {
            this->matrix[i] = unique_ptr<double[]>(new double[this->m]);
            for (int j = 0; j < matrix.m; j++) {
                this->matrix[i][j] = matrix.matrix[i][j];
            }
        }
    }
};


class SquareMatrix : public Matrix {
public:
    Matrix operator+(const Matrix& matrix) {
        Matrix first(*this);
        return first + matrix;
    }
    Matrix operator-(const Matrix& matrix) {
        Matrix first(*this);
        return first - matrix;
    }
    Matrix operator*(const Matrix& matrix) {
        Matrix first(*this);
        return first * matrix;
    }
    Matrix* operator+(Matrix* matrix) {
        return (SquareMatrix*)(*(Matrix*)this + matrix);
    }
    Matrix* operator-(Matrix* matrix) {
        return (SquareMatrix*)(*(Matrix*)this - matrix);
    }
    Matrix* operator*(Matrix* matrix) {
        return (SquareMatrix*)(*(Matrix*)this * matrix);
    }
    Matrix operator=(const Matrix& matrix) {
        if (matrix.getN() != matrix.getM()) {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
        this->n = matrix.getN();
        this->m = matrix.getM();
        this->matrix.release();
        this->matrix = unique_ptr<unique_ptr<double[]>[]>(new unique_ptr<double[]>[this->n]);
        for (int i = 0; i < matrix.getN(); i++) {
            this->matrix[i] = unique_ptr<double[]>(new double[this->m]);
            for (int j = 0; j < matrix.getM(); j++) {
                this->matrix[i][j] = matrix.getElem(i, j);
            }
        }
        return *this;
    }
    SquareMatrix() = delete;
    SquareMatrix(int n) : Matrix(n, n) {};
    SquareMatrix(const Matrix& matrix) : Matrix(matrix) {
        if (matrix.getN() != matrix.getM()) {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
    }
};


class IdentityMatrix : public SquareMatrix {
public:
    IdentityMatrix() = delete;
    IdentityMatrix(int n) : SquareMatrix(n) {
        for (int i = 0; i < n; i++) {
            this->matrix[i][i] = 1;
        }
    }
};


class PermutationMatrix : public SquareMatrix {
public:
    PermutationMatrix() = delete;
    PermutationMatrix(const Matrix& matrix, int v, int u) : SquareMatrix(matrix.getN()) {
        if (v > this->n or u > this->m) {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
        for (int i = 0; i < this->n; i++) {
            if (i == v - 1) {
                this->matrix[i][u - 1] = 1;
            }
            else if (i == u - 1) {
                this->matrix[i][v - 1] = 1;
            }
            else {
                this->matrix[i][i] = 1;
            }
        }
    }
};


class EliminationMatrix : public SquareMatrix {
public:
    EliminationMatrix() = delete;
    EliminationMatrix(const Matrix& matrix, int e_v, int e_u, int p_v, int p_u) : SquareMatrix(matrix.getN()) {
        if (e_v > this->n or e_u > matrix.getM() or p_v > this->n or p_u > matrix.getM()) {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
        for (int i = 0; i < this->n; i++) {
            this->matrix[i][i] = 1;
            if (i == e_v - 1) {
                this->matrix[i][p_v - 1] = -matrix.getElem(i, e_u - 1) / matrix.getElem(p_v - 1, p_u - 1);
            }
        }
    }
};


class NormalizationMatrix : public SquareMatrix {
public:
    NormalizationMatrix() = delete;
    NormalizationMatrix(const Matrix& matrix, int row) : SquareMatrix(matrix.getN()) {
        for (int i = 0; i < n; i++) {
            this->matrix[i][i] = 1;
            if (row - 1 == i) {
                this->matrix[i][i] = 1 / matrix.getElem(i, i);
            }
        }
    }
};


class ColumnVector : public Matrix {
public:
    double norm() {
        double res = 0;
        for (int i = 0; i < this->n; i++) {
            res += pow(this->matrix[i][0], 2);
        }
        return sqrt(res);
    }
    ColumnVector() = delete;
    ColumnVector(int n) : Matrix(n, 1) { }
    ColumnVector(int n, double val) : Matrix(n, 1) {
        for (int i = 0; i < n; i++) {
            this->matrix[i][0] = val;
        }
    }
    ColumnVector(const Matrix& matrix) : Matrix(matrix.getN(), 1) {
        if (matrix.getM() != 1) {
            throw invalid_argument("Error: the dimensional problem occurred\n");
        }
        for (int i = 0; i < matrix.getN(); i++) {
            this->matrix[i][0] = matrix.getElem(i, 0);
        }
    }
    ColumnVector(const vector<double>& input) : Matrix(input.size(), 1) {
        for (int i = 0; i < input.size(); i++) {
            this->matrix[i][0] = input[i];
        }
    }
};


class DiagonalMatrix : public Matrix {
public:
    DiagonalMatrix() = delete;
    DiagonalMatrix(const ColumnVector& init) : Matrix(init.getN(), init.getN()) {
        for (int i = 0; i < this->n; i++) {
            this->matrix[i][i] = init.getElem(i, 0);
        }
    }
};


// Function that makes an upper triangular matrix from the given one
SquareMatrix makeUpperTrinagular(Matrix& matrix, Matrix& e, int& steps) {
    int pivots = 0;
    for (int j = 0; j < matrix.getM(); j++) {
        if (pivots >= matrix.getN() - 1) {
            break;
        }
        double maxVal = abs(matrix.getElem(pivots, j));
        int row = pivots;
        for (int i = pivots + 1; i < matrix.getN(); i++) {
            if (abs(matrix.getElem(i, j)) > maxVal) {
                maxVal = abs(matrix.getElem(i, j));
                row = i;
            }
        }
        if (maxVal == 0) {
            continue;
        }
        if (row != pivots) {
            PermutationMatrix tmpMatrix(matrix, pivots + 1, row + 1);
            matrix = tmpMatrix * matrix;
            e = tmpMatrix * e;
        }
        for (int i = pivots + 1; i < matrix.getN(); i++) {
            if (matrix.getElem(i, j) != 0) {
                EliminationMatrix tmpMatrix(matrix, i + 1, j + 1, pivots + 1, j + 1);
                matrix = tmpMatrix * matrix;
                e = tmpMatrix * e;
            }
        }
        pivots++;
    }
    return e;
}


// Function that makes a down triangular matrix from the given one
SquareMatrix makeDownTrinagular(Matrix& matrix, Matrix& e, int& steps, bool second) {
    int pivots = 0;
    for (int j = matrix.getM() - 1; j > 0; j--) {
        if (pivots >= matrix.getN() - 1) {
            break;
        }
        if (!second) {
            double maxVal = abs(matrix.getElem(matrix.getN() - pivots - 1, j));
            int row = matrix.getN() - pivots - 1;
            for (int i = matrix.getN() - pivots - 2; i >= 0; i--) {
                if (abs(matrix.getElem(i, j)) > maxVal) {
                    maxVal = abs(matrix.getElem(i, j));
                    row = i;
                }
            }
            if (maxVal == 0) {
                continue;
            }
            if (row != matrix.getN() - pivots - 1) {
                PermutationMatrix tmpMatrix(matrix, matrix.getN() - pivots, row + 1);
                matrix = tmpMatrix * matrix;
                e = tmpMatrix * e;
            }
        }
        for (int i = matrix.getN() - pivots - 2; i >= 0; i--) {
            if (matrix.getElem(i, j) != 0) {
                EliminationMatrix tmpMatrix(matrix, i + 1, j + 1, matrix.getN() - pivots, j + 1);
                matrix = tmpMatrix * matrix;
                e = tmpMatrix * e;
            }
        }
        pivots++;
    }
    return e;
}


// Function that checks if the matrix is singular
bool isSingular(SquareMatrix matrix) {
    int c = 0;
    SquareMatrix tmp(IdentityMatrix(matrix.getN()));
    makeUpperTrinagular(matrix, tmp, c);
    for (int i = 0; i < matrix.getN(); i++) {
        if (matrix.getElem(i, i) == 0) {
            return true;
        }
    }
    return false;
}


// Function that normalizes a diagonal matrix
SquareMatrix diagonalNormalization(SquareMatrix& matrix, Matrix& e) {
    for (int i = 0; i < matrix.getN(); i++) {
        if (matrix.getElem(i, i) == 0) {
            continue;
        }
        NormalizationMatrix tmpMatrix(matrix, i + 1);
        matrix = tmpMatrix * matrix;
        e = tmpMatrix * e;
    }
    return e;
}


// Function that calculates inverse matrix
SquareMatrix findInverse(SquareMatrix matrix) {
    int c = 0;
    if (isSingular(matrix)) {
        throw invalid_argument("Error: matrix A is singular\n");
    }
    SquareMatrix e(IdentityMatrix(matrix.getN()));
    makeUpperTrinagular(matrix, e, c);
    makeDownTrinagular(matrix, e, c, true);
    diagonalNormalization(matrix, e);
    return e;
}


// Function that checks if the give transportation problem is balanced
void validityCheck(const Matrix& s, const Matrix& d) {
    double balance = 0;
    for (int i = 0; i < s.getN(); i++) {
        balance += s.getElem(i, 0);
    }
    for (int j = 0; j < d.getM(); j++) {
        balance -= d.getElem(0, j);
    }
    if (balance != 0) {
        throw invalid_argument("Error: The problem is not balanced!\n");
    }
}


// Function that implements North-West method
double northWest(Matrix s, Matrix c, Matrix d) {
    int i = 0, j = 0;
    double result = 0;
    cout << "North - West corner method approximation:\nBasic variables vector: ";

    // Until we do not reach the right-bottom corner
    while (i < c.getN() && j < c.getM()) {
        cout << "(" << i << ", " << j << "); ";

        // If the demand is higher or equal to the supply
        if (d.getElem(0, j) >= s.getElem(i, 0)) {
            result += c.getElem(i, j) * s.getElem(i, 0);
            d.setElem(0, j, d.getElem(0, j) - s.getElem(i, 0));
            i++;
        }
        // If the supply is higher than the demand
        else {
            result += c.getElem(i, j) * d.getElem(0, j);
            s.setElem(i, 0, s.getElem(i, 0) - d.getElem(0, j));
            j++;
        }
    }
    cout << "\nSolution: " << result << "\n";
    return result;
}

// Function that implements Vogel's method
double vogel(Matrix S, Matrix C, Matrix D) {
    int n = C.getN(), m = C.getM();
    vector<double> rowDifferencies(n), columnDifferencies(m);
    Matrix ans(n, m);
    double result = 0;
    cout << "\nVogel's method approximation: \nInitial solution vector: \n";
    while (true) {
        // Calculating row differencies
        for (int i = 0; i < n; i++) {
            double minim = 9999, secMinim = 9999;
            int minimIndex = 9999;
            double curr_value;
            // Finding minimum
            for (int j = 0; j < m; j++) {
                curr_value = C.getElem(i, j);
                if (curr_value < minim && curr_value != -1) {
                    minim = curr_value;
                    minimIndex = j;
                }
            }

            // Finding second minimum
            for (int j = 0; j < m; j++) {
                curr_value = C.getElem(i, j);
                if (curr_value < secMinim && curr_value >= minim && minimIndex != j && curr_value != -1) {
                    secMinim = curr_value;
                }
            }
            if (secMinim == 9999 && minim == 9999) rowDifferencies[i] = -1;
            else if (secMinim == 9999) rowDifferencies[i] = minim;
            else rowDifferencies[i] = secMinim - minim;
        }

        // Calculating column differencies
        for (int i = 0; i < m; i++) {
            double minim = 9999, secMinim = 9999;
            int minimIndex = 9999;
            double curr_value;
            for (int j = 0; j < n; j++) {
                curr_value = C.getElem(j, i);
                if (curr_value < minim && curr_value != -1) {
                    minim = curr_value;
                    minimIndex = j;
                }
            }
            for (int j = 0; j < n; j++) {
                curr_value = C.getElem(j, i);
                if (curr_value < secMinim && curr_value >= minim && minimIndex != j && curr_value != -1) {
                    secMinim = curr_value;
                }
            }
            if (secMinim == 9999 || minim == 9999) columnDifferencies[i] = -1;
            else if (secMinim == 9999) columnDifferencies[i] = minim;
            else columnDifferencies[i] = secMinim - minim;
        }
        // Defining maximum value
        double maxRow = -1, maxColumn = -1;
        int initialRowIndex = -1, initialColumnIndex = -1;
        for (int i = 0; i < n; i++) {
            double value = rowDifferencies[i];
            if (value > maxRow) {
                maxRow = value;
                initialRowIndex = i;
            }
        }
        for (int i = 0; i < m; i++) {
            double value = columnDifferencies[i];
            if (value > maxColumn) {
                maxColumn = value;
                initialColumnIndex = i;
            }
        }

        if (max(maxRow, maxColumn) == -1) {
            break;
        }

        int rowIndex, columnIndex;
        double minValue = 99999999;

        // If row is maximum
        if (max(maxRow, maxColumn) == maxRow) {
            rowIndex = initialRowIndex;
            columnIndex;
            // Defining min cell
            for (int i = 0; i < m; i++) {
                if (C.getElem(rowIndex, i) < minValue && C.getElem(rowIndex, i) != -1) {
                    minValue = C.getElem(rowIndex, i);
                    columnIndex = i;
                }
            }
        }
        // if column is maximum
        else {
            rowIndex = -1;
            columnIndex = initialColumnIndex;
            // Defining min cell
            for (int i = 0; i < n; i++) {
                if (C.getElem(i, columnIndex) < minValue && C.getElem(i, columnIndex) != -1) {
                    minValue = C.getElem(i, columnIndex);
                    rowIndex = i;
                }
            }
            
        }
        // If demand is over delete column
        if (D.getElem(0, columnIndex) <= S.getElem(rowIndex, 0)) {
            ans.setElem(rowIndex, columnIndex, D.getElem(0, columnIndex));
            result += C.getElem(rowIndex, columnIndex) * D.getElem(0, columnIndex);
            S.setElem(rowIndex, 0, S.getElem(rowIndex, 0) - D.getElem(0, columnIndex));
            for (int i = 0; i < n; i++) {
                C.setElem(i, columnIndex, -1);
            }
        }
        // if supply is over delete row
        else {
            ans.setElem(rowIndex, columnIndex, S.getElem(rowIndex, 0));
            result += C.getElem(rowIndex, columnIndex) * S.getElem(rowIndex, 0);
            D.setElem(0, columnIndex, D.getElem(0, columnIndex) - S.getElem(rowIndex, 0));
            for (int i = 0; i < m; i++) {
                C.setElem(rowIndex, i, -1);
            }
        }
    }
    cout << ans;
    cout << "Total distribution cost: " << result;
    return result;
}

// Function that implements Russel's method
double russel(Matrix s, Matrix c, Matrix d) {
    vector<double> u(s.getN());
    vector<double> v(d.getM());
    double result = 0;
    cout << "Russel's method approximation:\nBasic variables vector:";

    // The cycle stops when there are no rows or columns left
    while (true) {
        fill(u.begin(), u.end(), -1);
        fill(v.begin(), v.end(), -1);

        // Finding the largest values in rows and columns for vetors U and V
        for (int i = 0; i < s.getN(); i++) {
            for (int j = 0; j < d.getM(); j++) {
                u[i] = max(u[i], c.getElem(i, j));
                v[j] = max(v[j], c.getElem(i, j));
            }
        }

        // Finding the proper delta to consider
        double delta = -1;
        int deltaI = -1, deltaJ = -1;
        for (int i = 0; i < s.getN(); i++) {
            for (int j = 0; j < d.getM(); j++) {
                if (u[i] < 0 || v[i] < 0) {
                    continue;
                }
                double tmp = c.getElem(i, j) - (u[i] + v[j]);
                if (tmp < 0 && abs(tmp) > delta) {
                    delta = abs(tmp);
                    deltaI = i;
                    deltaJ = j;
                }
            }
        }

        // If no delta found, then we are out of rows and columns
        if (delta == -1) {
            break;
        }
        cout << "(" << deltaI << ", " << deltaJ << "); ";

        // If the demand is higher or equal to the supply 
        if (d.getElem(0, deltaJ) >= s.getElem(deltaI, 0)) {
            result += c.getElem(deltaI, deltaJ) * s.getElem(deltaI, 0);
            for (int j = 0; j < d.getM(); j++) {
                c.setElem(deltaI, j, -1);
            }
        }
        // If the supply is higher than the demand
        else {
            result += c.getElem(deltaI, deltaJ) * d.getElem(0, deltaJ);
            for (int i = 0; i < s.getN(); i++) {
                c.setElem(i, deltaJ, -1);
            }
        }
    }
    cout << "\nSolution: " << result << "\n";
    return result;
}


int main() {
    try {
        int n, m;
        cin >> n >> m;
        Matrix s(1, n), c(n, m), d(1, m);
        cin >> s >> c >> d;
        s = s.transpose();
        validityCheck(s, d);
        northWest(s, c, d);
        vogel(s, c, d);
        russel(s, c, d);
    }
    catch (const exception& ex) {
        cout << ex.what() << "\n";
    }
}


/*
Input example:


3 5
140 180 160
2 3 4 2 4
8 4 1 4 1
9 7 3 7 2
60 70 120 130 100
*/
