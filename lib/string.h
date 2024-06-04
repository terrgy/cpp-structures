#include <cstring>
#include <algorithm>
#include <iostream>
#include <cctype>

class String {
private:
    size_t p_size;
    size_t p_capacity;
    char* p_data;

    void _change_capacity(size_t new_capacity) {
        char* new_data = new char[new_capacity + 1];
        std::copy(p_data, p_data + p_size + 1, new_data);
        delete[] p_data;
        p_data = new_data;
        p_capacity = new_capacity;
    }

    bool _check_substr(size_t start, const String& substr) const {
        for (size_t position = 0; position < substr.p_size; ++position) {
            if (p_data[start + position] != substr.p_data[position]) {
                return false;
            }
        }
        return true;
    }

    size_t _find(const String& substr, size_t search_start, size_t search_end, int modifier) const {
        if (p_size < substr.p_size) {
            return length();
        }

        for (size_t start = search_start; ; start = static_cast<size_t>(static_cast<int>(start) + modifier)) {
            if (_check_substr(start, substr)) {
                return start;
            }
            if (start == search_end) {
                return length();
            }
        }
    }

    void _set_last_null_char() {
        p_data[p_size] = '\0';
    }

    String(size_t size, char* data) : p_size(size), p_capacity(size), p_data(data) {}

    void _swap(String& other) noexcept {
        std::swap(p_size, other.p_size);
        std::swap(p_capacity, other.p_capacity);
        std::swap(p_data, other.p_data);
    }

public:
    String() : p_size(0), p_capacity(0), p_data(new char[1]) {
        p_data[0] = '\0';
    }

    String(const char* str) :
            p_size(strlen(str)),
            p_capacity(p_size),
            p_data(new char[p_size + 1])
    {
        std::copy(str, str + p_size, p_data);
        _set_last_null_char();
    }

    String(size_t size, char fill_char) : p_size(size), p_capacity(size), p_data(new char[size + 1]) {
        std::fill(p_data, p_data + p_size, fill_char);
        _set_last_null_char();
    }

    String(const String& other) : p_size(other.p_size), p_capacity(p_size), p_data(new char[p_size + 1]) {
        std::copy(other.p_data, other.p_data + p_size + 1, p_data);
    }

    String& operator=(const String& other) {
        if (this != &other) {
            if (p_capacity >= other.p_size) {
                std::copy(other.p_data, other.p_data + other.p_size + 1, p_data);
                p_size = other.p_size;
            } else {
                String tmp(other);
                _swap(tmp);
            }
        }
        return *this;
    }

    ~String() {
        delete[] p_data;
    }

    size_t size() const {
        return p_size;
    }

    size_t capacity() const {
        return p_capacity;
    }

    size_t length() const {
        return size();
    }

    char* data() {
        return p_data;
    }

    const char* data() const {
        return p_data;
    }

    char& operator[](size_t index) {
        return p_data[index];
    }

    const char& operator[](size_t index) const {
        return p_data[index];
    }

    void pop_back() {
        --p_size;
        _set_last_null_char();
    }

    void push_back(char new_char) {
        if (p_size == p_capacity) {
            _change_capacity(p_capacity * 2 + 1);
        }

        p_data[p_size++] = new_char;
        _set_last_null_char();
    }

    char& front() {
        return p_data[0];
    }

    const char& front() const {
        return p_data[0];
    }

    char& back() {
        return p_data[p_size - 1];
    }

    const char& back() const {
        return p_data[p_size - 1];
    }

    String& operator+=(char new_char) {
        push_back(new_char);
        return *this;
    }

    String& operator+=(const String& other) {
        if (p_size + other.p_size > p_capacity) {
            _change_capacity(p_size + other.p_size);
        }

        std::copy(other.p_data, other.p_data + other.p_size + 1, p_data + p_size);
        p_size += other.p_size;
        return *this;
    }

    void shrink_to_fit() {
        if (p_size == p_capacity) {
            return;
        }

        _change_capacity(p_size);
    }

    bool empty() const {
        return p_size == 0;
    }

    void clear() {
        p_size = 0;
        _set_last_null_char();
    }

    String substr(size_t start, size_t count) const {
        if (p_size - start < count) {
            count = p_size - start;
        }
        char* new_data = new char[count + 1];
        std::copy(p_data + start, p_data + start + count, new_data);
        new_data[count] = '\0';
        return String(count, new_data);
    }

    size_t find(const String& substr) const {
        return _find(substr, 0, p_size - substr.p_size, 1);
    }

    size_t rfind(const String& substr) const {
        return _find(substr, p_size - substr.p_size, 0, -1);
    }

    void reserve(size_t capacity) {
        if (p_capacity < capacity) {
            _change_capacity(capacity);
        }
    }
};

bool operator==(const String& left_str, const String& right_str) {
    return (left_str.size() == right_str.size()) && (strcmp(left_str.data(), right_str.data()) == 0);
}

bool operator!=(const String& left_str, const String& right_str) {
    return !(left_str == right_str);
}

bool operator<(const String& left_str, const String& right_str) {
    return strcmp(left_str.data(), right_str.data()) < 0;
}

bool operator>(const String& left_str, const String& right_str) {
    return right_str < left_str;
}

bool operator<=(const String& left_str, const String& right_str) {
    return !(right_str < left_str);
}

bool operator>=(const String& left_str, const String& right_str) {
    return !(left_str < right_str);
}

String operator+(const String& left_str, const String& right_str) {
    String result;
    result.reserve(left_str.size() + right_str.size());
    result += left_str;
    result += right_str;
    return result;
}

String operator+(const String& str, char new_char) {
    String result;
    result.reserve(str.size() + 1);
    result += str;
    result.push_back(new_char);
    return result;
}

String operator+(char left_value, const String& right_value) {
    String result;
    result.reserve(right_value.size() + 1);
    result.push_back(left_value);
    result += right_value;
    return result;
}

std::istream& operator>>(std::istream& input, String& str) {
    str.clear();
    char new_char;
    while (input.get(new_char)) {
        if (std::isspace(new_char) == 0) {
            str.push_back(new_char);
            break;
        }
    }
    while (input.get(new_char)) {
        if (std::isspace(new_char) != 0) {
            input.unget();
            break;
        }
        str.push_back(new_char);
    }
    return input;
}

std::ostream& operator<<(std::ostream& output, const String& str) {
    output << str.data();
    return output;
}