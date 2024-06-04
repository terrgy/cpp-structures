#include <utility>
#include <vector>
#include <string>
#include <cassert>
#include <compare>
#include <complex>
#include <cmath>
#include <algorithm>

class BigInteger {
public:
    enum class Sign {
        Negative = -1,
        Zero = 0,
        Positive = 1
    };

    static std::strong_ordering absCompare(const BigInteger& number1, const BigInteger& number2) {
        auto size_ordering = number1.digits.size() <=> number2.digits.size();
        if (size_ordering != std::strong_ordering::equal) return size_ordering;
        for (size_t i = number1.digits.size() - 1;; --i) {
            auto digit_ordering = number1.digits[i] <=> number2.digits[i];
            if (digit_ordering != std::strong_ordering::equal) return digit_ordering;
            if (i == 0) {
                break;
            }
        }
        return std::strong_ordering::equal;
    }

private:
    using fft_type = std::complex<double>;
    static const int DECIMAL_BASE = 1000;
    static const size_t DECIMAL_BASE_COUNT = 3;
    std::vector<int> digits = {0};
    Sign sign = Sign::Zero;

    void _normalize() {
        while ((digits.size() > 1) && (digits.back() == 0)) {
            digits.pop_back();
        }
        if ((digits.size() == 1) && (digits[0] == 0)) {
            sign = Sign::Zero;
        }
    }

    static Sign _multiplySigns(Sign sign1, Sign sign2) {
        return static_cast<Sign>(static_cast<int>(sign1) * static_cast<int>(sign2));
    }

    static std::pair<int, int> _digitDivMod(int a, int multiplier) {
        int div = a / DECIMAL_BASE, mod = a % DECIMAL_BASE;
        if (mod * multiplier < 0) {
            div -= multiplier;
            mod += DECIMAL_BASE * multiplier;
        }
        return {div, mod};
    }

    Sign _determineAddSign(const BigInteger& other, Sign other_sign) {
        std::strong_ordering compare_res = absCompare(*this, other);

        if (compare_res != std::strong_ordering::equal) {
            return compare_res == std::strong_ordering::less ? other_sign : sign;
        }
        return sign == other_sign ? sign : Sign::Zero;
    }

    void _selfAddition(const BigInteger& other, int multiplier = 1) {
        Sign other_sign = static_cast<Sign>(static_cast<int>(other.sign) * multiplier);
        Sign result_sign = _determineAddSign(other, other_sign);
        if (result_sign == Sign::Zero) {
            makeZero();
            return;
        }

        digits.resize(std::max(digits.size(), other.digits.size()));
        int carry = 0;
        for (size_t i = 0; (i < other.digits.size()) || (carry && (i < digits.size())); ++i) {
            int cur_digit_sum = static_cast<int>(sign) * digits[i] + carry;
            if (i < other.digits.size()) {
                cur_digit_sum += static_cast<int>(other_sign) * other.digits[i];
            }
            auto [new_carry, new_digit] = _digitDivMod(cur_digit_sum, static_cast<int>(result_sign));
            carry = new_carry;
            digits[i] = abs(new_digit);
        }
        if (carry) {
            digits.push_back(abs(carry));
        }

        sign = result_sign;
        _normalize();
    }

    static void _fftReverseBits(std::vector<fft_type>& f) {
        for (size_t i = 1, j = 0; i < f.size(); ++i) {
            size_t bit = f.size() >> 1;
            for (; j >= bit; bit >>= 1) {
                j -= bit;
            }
            j += bit;
            if (i < j) {
                swap(f[i], f[j]);
            }
        }
    }

    static void _fftBase(std::vector<fft_type>& f, bool invert = false) {
        _fftReverseBits(f);

        for (size_t len = 2; len <= f.size(); len <<= 1) {
            double ang = 2 * M_PI / static_cast<double>(len) * (invert ? -1 : 1);
            fft_type w_move(cos(ang), sin(ang));
            for (size_t i = 0; i < f.size(); i += len) {
                fft_type w(1);
                for (size_t j = 0; j < len / 2; ++j) {
                    fft_type u = f[i + j], v = f[i + j + len / 2] * w;
                    f[i + j] = u + v;
                    f[i + j + len / 2] = u - v;
                    w *= w_move;
                }
            }
        }
    }

    static void _fft(std::vector<fft_type>& f) {
        _fftBase(f);
    }

    static void _invFft(std::vector<fft_type>& f) {
        _fftBase(f, true);
        for (size_t i = 0; i < f.size(); ++i)
        {
            f[i] /= static_cast<double>(f.size());
        }
    }

    static int _findQuotDigit(const BigInteger& remainder, const BigInteger& divider) {
        int bin_search_left = 0, bin_search_right = DECIMAL_BASE;
        while (bin_search_right - bin_search_left > 1) {
            int bin_search_mid = (bin_search_left + bin_search_right) / 2;
            BigInteger sub = divider;
            sub *= bin_search_mid;
            if (absCompare(sub, remainder) == std::strong_ordering::greater) {
                bin_search_right = bin_search_mid;
            } else {
                bin_search_left = bin_search_mid;
            }
        }
        return bin_search_left;
    }

    static std::pair<BigInteger, BigInteger> _bigintDivMod(const BigInteger& number1, const BigInteger& number2) {
        BigInteger remainder, quot;
        if ((number1.digits.size() < number2.digits.size()) || number1.isZero()) {
            quot.makeZero();
            remainder = number1;
            return {quot, remainder};
        }

        remainder.sign = Sign::Positive;
        remainder.digits.resize(number2.digits.size());
        std::copy(number1.digits.end() - static_cast<int>(number2.digits.size()), number1.digits.end(),
                  remainder.digits.begin());
        std::vector<int> quot_digits;
        quot_digits.reserve(number1.digits.size() - number2.digits.size() + 1);

        for (int i = static_cast<int>(number1.digits.size() - number2.digits.size() - 1); i >= -1; --i) {
            int quot_digit = _findQuotDigit(remainder, number2);
            quot_digits.push_back(quot_digit);

            BigInteger sub = number2;
            sub.sign = Sign::Positive;
            remainder -= (sub *= quot_digit);
            if (i != -1) {
                remainder.digits.insert(remainder.digits.begin(),number1.digits[static_cast<size_t>(i)]);
                remainder.sign = Sign::Positive;
                remainder._normalize();
            }
        }

        std::reverse(quot_digits.begin(), quot_digits.end());
        quot.digits = std::move(quot_digits);
        quot.sign = _multiplySigns(number1.sign, number2.sign);
        quot._normalize();

        remainder.sign = number1.sign;
        remainder._normalize();
        return {quot, remainder};
    }

    void _multiplySmall(int number) {
        int carry = 0;
        for (int & digit : digits) {
            digit = digit * number + carry;
            carry = digit / DECIMAL_BASE;
            digit %= DECIMAL_BASE;
        }
        if (carry) {
            digits.push_back(carry);
        }
    }

public:
    BigInteger() = default;

    explicit BigInteger(const std::string& str) : sign(Sign::Positive) {
        size_t negative_pos_mod = 0;
        if (str.front() == '-') {
            sign = Sign::Negative;
            negative_pos_mod = 1;
        }
        digits.resize((str.size() - 1 - negative_pos_mod) / DECIMAL_BASE_COUNT + 1);
        for (size_t i = 0; i < digits.size() - 1; ++i) {
            digits[i] = std::stoi(str.substr(str.size() - (i + 1) * DECIMAL_BASE_COUNT, DECIMAL_BASE_COUNT));
        }
        digits.back() = std::stoi(str.substr(negative_pos_mod,
                                             (str.size() - 1 - negative_pos_mod) % DECIMAL_BASE_COUNT + 1));

        _normalize();
    }

    BigInteger(long long value) {
        if (!value) {
            return;
        }
        if (value < 0) {
            sign = Sign::Negative;
            value *= -1;
        } else {
            sign = Sign::Positive;
        }
        digits = {};
        while (value) {
            digits.push_back(static_cast<int>(value % DECIMAL_BASE));
            value /= DECIMAL_BASE;
        }
    }

    std::string toString() const {
        std::string result;
        result.reserve(digits.size() * DECIMAL_BASE_COUNT + 1);
        if (isNegative()) {
            result += '-';
        }
        result += std::to_string(digits.back());
        for (int i = static_cast<int>(digits.size()) - 2; i > -1; --i) {
            std::string part_digit = std::to_string(digits[static_cast<size_t>(i)]);
            result.append(DECIMAL_BASE_COUNT - part_digit.size(), '0');
            result += part_digit;
        }
        return result;
    }

    explicit operator bool() const {
        return !isZero();
    }

    BigInteger& operator+=(const BigInteger& other) {
        _selfAddition(other);
        return *this;
    }

    BigInteger& operator-=(const BigInteger& other) {
        _selfAddition(other, -1);
        return *this;
    }

    void changeSign() {
        sign = static_cast<Sign>(-static_cast<int>(sign));
    }

    BigInteger operator-() const {
        BigInteger result(*this);
        result.changeSign();
        return result;
    }

    BigInteger& operator++() {
        (*this) += 1;
        return *this;
    }

    BigInteger operator++(int) {
        BigInteger result = *this;
        ++(*this);
        return result;
    }

    BigInteger& operator--() {
        (*this) -= 1;
        return *this;
    }

    BigInteger operator--(int) {
        BigInteger result = *this;
        --(*this);
        return result;
    }

    BigInteger& operator*=(const BigInteger& other) {
        sign = _multiplySigns(sign, other.sign);
        if (isZero()) {
            makeZero();
            return *this;
        }

        if (other.digits.size() == 1) {
            _multiplySmall(other.digits[0]);
            return *this;
        }

        size_t target_size = digits.size() + other.digits.size();
        size_t f_size = 1;
        while (f_size < target_size)
        {
            f_size *= 2;
        }

        std::vector<fft_type> number1_f(f_size), number2_f(f_size);
        copy(digits.begin(), digits.end(), number1_f.begin());
        copy(other.digits.begin(), other.digits.end(), number2_f.begin());

        _fft(number1_f);
        _fft(number2_f);

        for (size_t i = 0; i < f_size; ++i)
        {
            number1_f[i] *= number2_f[i];
        }

        _invFft(number1_f);

        digits.resize(f_size);
        long long carry = 0;
        for (size_t i = 0; i < f_size; ++i)
        {
            long long digit = (carry + static_cast<long long>(round(number1_f[i].real())));
            carry = digit / DECIMAL_BASE;
            digits[i] = static_cast<int>(digit % DECIMAL_BASE);
        }
        _normalize();
        return *this;
    }

    BigInteger& operator/=(const BigInteger& other) {
        (*this) = _bigintDivMod(*this, other).first;
        return *this;
    }

    BigInteger& operator%=(const BigInteger& other) {
        (*this) = _bigintDivMod(*this, other).second;
        return *this;
    }

    Sign getSign() const {
        return sign;
    }

    void multiplyDecimals(size_t decimals) {
        static const int ONE_DECIMAL_DIGIT = 10;

        if (isZero()) {
            return;
        }
        if (decimals % DECIMAL_BASE_COUNT) {
            size_t move_decimals_count = decimals % DECIMAL_BASE_COUNT;
            int decimals_multiplier = 1;
            for (size_t i = 0; i < move_decimals_count; ++i) {
                decimals_multiplier *= ONE_DECIMAL_DIGIT;
            }
            int carry = 0;
            for (int & digit : digits) {
                long long cur_digit = digit * decimals_multiplier + carry;
                carry = static_cast<int>(cur_digit / DECIMAL_BASE);
                digit = static_cast<int>(cur_digit % DECIMAL_BASE);
            }
            if (carry) {
                digits.push_back(carry);
            }
        }
        digits.insert(digits.begin(), decimals / DECIMAL_BASE_COUNT, 0);
    }

    bool isZero() const {
        return sign == Sign::Zero;
    }

    bool isPositive() const {
        return sign == Sign::Positive;
    }

    bool isNegative() const {
        return sign == Sign::Negative;
    }

    void makeZero() {
        digits = {0};
        sign = Sign::Zero;
    }

    void makePositive() {
        if (!isZero()) {
            sign = Sign::Positive;
        }
    }

    void makeNegative() {
        if (!isZero()) {
            sign = Sign::Negative;
        }
    }

    bool operator==(const BigInteger&) const = default;
};

std::strong_ordering operator<=>(const BigInteger& number1, const BigInteger& number2) {
    auto sign_ordering = number1.getSign() <=> number2.getSign();
    if (sign_ordering != std::strong_ordering::equal) {
        return sign_ordering;
    }
    return number1.isNegative() ? BigInteger::absCompare(number2, number1) : BigInteger::absCompare(number1, number2);
}

BigInteger operator%(const BigInteger& number1, const BigInteger& number2) {
    BigInteger result = number1;
    result %= number2;
    return result;
}

BigInteger operator/(const BigInteger& number1, const BigInteger& number2) {
    BigInteger result = number1;
    result /= number2;
    return result;
}

BigInteger operator*(const BigInteger& number1, const BigInteger& number2) {
    BigInteger result = number1;
    result *= number2;
    return result;
}

BigInteger operator+(const BigInteger& number1, const BigInteger& number2) {
    BigInteger result = number1;
    result += number2;
    return result;
}

BigInteger operator-(const BigInteger& number1, const BigInteger& number2) {
    BigInteger result = number1;
    result -= number2;
    return result;
}

std::istream& operator>>(std::istream& in, BigInteger& number) {
    std::string str;
    in >> str;
    number = BigInteger(str);
    return in;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& number) {
    out << number.toString();
    return out;
}

BigInteger operator""_bi(unsigned long long number) {
    return {static_cast<long long>(number)};
}

BigInteger operator""_bi(const char* str) {
    return BigInteger(str);
}

class Rational {
private:
    BigInteger numerator = 0;
    BigInteger denominator = 1;

    void _normalizeSign() {
        if (denominator.isNegative()) {
            numerator.changeSign();
            denominator.changeSign();
        }
    }

    void _normalizePrime()  {
        BigInteger gcd1 = numerator;
        BigInteger gcd2 = denominator;
        gcd1.makePositive();
        while (gcd2) {
            gcd1 %= gcd2;
            std::swap(gcd1, gcd2);
        }
        numerator /= gcd1;
        denominator /= gcd1;
    }

public:
    Rational() = default;
    Rational(BigInteger number) : numerator(std::move(number)) {}
    Rational(int number) : numerator(number) {}
    Rational(BigInteger numerator, BigInteger denominator) :
            numerator(std::move(numerator)), denominator(std::move(denominator))
    {
        _normalizeSign();
        _normalizePrime();
    }

    Rational& operator+=(const Rational& other) {
        numerator = numerator * other.denominator + denominator * other.numerator;
        denominator *= other.denominator;
        _normalizePrime();
        return *this;
    }

    Rational& operator-=(const Rational& other) {
        numerator = numerator * other.denominator - denominator * other.numerator;
        denominator *= other.denominator;
        _normalizePrime();
        return *this;
    }

    Rational& operator*=(const Rational& other) {
        numerator *= other.numerator;
        denominator *= other.denominator;
        _normalizePrime();
        return *this;
    }

    Rational& operator/=(const Rational& other) {
        numerator *= other.denominator;
        denominator *= other.numerator;
        _normalizeSign();
        _normalizePrime();
        return *this;
    }

    void changeSign() {
        numerator.changeSign();
    }

    Rational operator-() const {
        Rational result = *this;
        result.changeSign();
        return result;
    }

    const BigInteger& getNumerator() const {
        return numerator;
    }

    const BigInteger& getDenominator() const {
        return denominator;
    }

    std::string toString() const {
        std::string result = numerator.toString();
        if (denominator != 1) {
            result += "/" + denominator.toString();
        }
        return result;
    }

    std::string asDecimal(size_t precision = 0) const {
        if (!precision) {
            return (numerator / denominator).toString();
        }
        BigInteger numerator_copy = numerator;
        numerator_copy.multiplyDecimals(precision);
        numerator_copy /= denominator;
        std::string result = numerator_copy.toString();
        bool is_negative = numerator_copy.isNegative();
        if (result.size() - is_negative <= precision) {
            result.insert(result.begin() + is_negative, precision + 1 - result.size() + is_negative, '0');
        }
        result.insert(result.end() - static_cast<int>(precision), '.');
        return result;
    }

    explicit operator double() const {
        static const int PRECISION = 20;
        return std::stod(asDecimal(PRECISION));
    }

    bool operator==(const Rational&) const = default;
};

std::strong_ordering operator<=>(const Rational& number1, const Rational& number2) {
    return (number1.getNumerator() * number2.getDenominator()) <=> (number1.getDenominator() * number2.getNumerator());
}

Rational operator+(const Rational& number1, const Rational& number2) {
    Rational result = number1;
    result += number2;
    return result;
}

Rational operator*(const Rational& number1, const Rational& number2) {
    Rational result = number1;
    result *= number2;
    return result;
}

Rational operator-(const Rational& number1, const Rational& number2) {
    Rational result = number1;
    result -= number2;
    return result;
}

Rational operator/(const Rational& number1, const Rational& number2) {
    Rational result = number1;
    result /= number2;
    return result;
}