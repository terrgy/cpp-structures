#include <array>
#include "biginteger.h"

template<size_t N, size_t D>
struct SqrtSearchHelper {
    static const int value = (D * D < N ? D : SqrtSearchHelper<N, D / 2>::value);
};

template<size_t N>
struct SqrtSearchHelper<N, 1> {
    static const int value = 1;
};

template<size_t N, size_t D>
struct PrimeSearchHelper {
    static const bool value = (N % D) && PrimeSearchHelper<N, D - 1>::value;
};

template<size_t N>
struct PrimeSearchHelper<N, 1> {
    static const bool value = true;
};

template<size_t N>
struct IsPrime {
    static const bool value = PrimeSearchHelper<N, 2 * SqrtSearchHelper<N, N - 1>::value >::value;
};

template<>
struct IsPrime<2> {
    static const bool value = true;
};

template<>
struct IsPrime<1> {
    static const bool value = false;
};

template<size_t N>
class Residue {
private:
    using ResidueType = unsigned long long;

    static_assert(N > 0);
    ResidueType value = 0;

    static ResidueType _binPow(ResidueType number, size_t power) {
        if (!power) {
            return 1;
        }
        ResidueType result = _binPow(number, power / 2);
        result = static_cast<ResidueType>(static_cast<unsigned long long>(result) * result % N);
        if (power % 2) {
            result = static_cast<ResidueType>(static_cast<unsigned long long>(result) * number % N);
        }
        return result;
    }

public:
    Residue() = default;
    explicit Residue(int value) {
        value %= N;
        if (value < 0) {
            value += N;
        }
        this->value = static_cast<ResidueType>(value);
    }

    Residue& operator+=(Residue another) {
        value += another.value;
        if (value >= N) {
            value -= N;
        }
        return *this;
    }

    Residue& operator-=(Residue another) {
        value -= another.value;
        if (value < 0) {
            value += N;
        }
        return *this;
    }

    Residue& operator*=(Residue another) {
        value = (value * another.value) % N;
        return *this;
    }

    explicit operator int() const {
        return static_cast<int>(value);
    }

    Residue& operator/=(Residue another) {
        static_assert(IsPrime<N>::value);
        *this *= _binPow(another, N - 2);
    }
};

template<size_t N>
Residue<N> operator+(Residue<N> left, Residue<N> right) {
    left += right;
    return left;
}

template<size_t N>
Residue<N> operator-(Residue<N> left, Residue<N> right) {
    left -= right;
    return left;
}

template<size_t N>
Residue<N> operator*(Residue<N> left, Residue<N> right) {
    left *= right;
    return left;
}

template<size_t N>
std::istream& operator>>(std::istream& in, Residue<N>& number) {
    int value;
    in >> value;
    number = Residue<N>(value);
    return in;
}

template<size_t N>
std::ostream& operator<<(std::ostream& out, Residue<N> number) {
    out << (int)number;
    return out;
}


template<size_t M, size_t N, typename Field=Rational>
class Matrix {
private:
    std::array< std::array<Field, N>, M> matrix;

public:
    Matrix() = default;
    explicit Matrix(const std::array< std:: array<Field, N>, M>& matrix) : matrix(matrix) {}

    Matrix& operator+=(const Matrix& another) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                matrix[i][j] += another.matrix[i][j];
            }
        }
        return *this;
    }

    Matrix& operator-=(const Matrix& another) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                matrix[i][j] -= another.matrix[i][j];
            }
        }
        return *this;
    }

    Matrix& operator*=(const Field& multiplier) {
        for (size_t i = 0; i < M; ++i) {
            for (size_t j = 0; j < N; ++j) {
                matrix[i][j] *= multiplier;
            }
        }
        return *this;
    }

    std::array<Field, N> getRow(size_t index) const {
        return matrix[index];
    }

    std::array<Field, M> getColumn(size_t index) const {
        std::array<Field, M> column;
        for (size_t i = 0; i < M; ++i) {
            column[i] = matrix[i][index];
        }
        return column;
    }

    void setRow(std::array<Field, N>& row, size_t index) {
        matrix[index] = row;
    }

    void setColumn(std::array<Field, M>& column, size_t index) {
        for (size_t i = 0; i < M; ++i) {
            matrix[i][index] = column[i];
        }
    }

    Matrix<N, M> transposed() const {
        Matrix<N, M> result;
        for (size_t i = 0; i < M; ++i) {
            result.setColumn(getColumn(i), i);
        }
    }

    Field trace() const {
        static_assert(N == M);
        Field result = 0;
        for (size_t i = 0; i < N; ++i) {
            result += matrix[i][i];
        }
        return result;
    }
};

template<size_t M, size_t N>
Matrix<M, N> operator+(const Matrix<M, N>& matrix1, const Matrix<M, N>& matrix2) {
    Matrix<M, N> result = matrix1;
    result += matrix2;
    return result;
}

template<size_t M, size_t N>
Matrix<M, N> operator-(const Matrix<M, N>& matrix1, const Matrix<M, N>& matrix2) {
    Matrix<M, N> result = matrix1;
    result -= matrix2;
    return result;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Matrix<M, N, Field>& matrix, const Field& multiplier) {
    Matrix<M, N, Field> result = matrix;
    result *= multiplier;
    return result;
}

template<size_t M, size_t N, typename Field>
Matrix<M, N, Field> operator*(const Field& multiplier, const Matrix<M, N, Field>& matrix) {
    return matrix * multiplier;
}
