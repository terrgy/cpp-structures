#include <iostream>
#include <memory>

template <size_t N>
class StackStorage {
private:
    char buffer[N];
    char* used;

public:
    StackStorage() : buffer(), used(buffer) {}

    char* allocate(size_t bytes, size_t alignment) {
        uintptr_t shift = reinterpret_cast<uintptr_t>(used) % alignment;
        char* new_used = used;
        if (shift) {
            new_used += alignment - shift;
        }
        if (new_used + bytes > buffer + N) {
            throw std::bad_alloc();
        }
        char* ptr = new_used;
        used = new_used + bytes;
        return ptr;
    }
};

template<typename T, size_t N>
class StackAllocator {
public:
    StackStorage<N>* storage;

    using value_type = T;

    explicit StackAllocator(const StackStorage<N>& storage) : storage( const_cast<StackStorage<N>*>(&storage) ) {}

    template<typename U>
    StackAllocator(const StackAllocator<U, N>& other) : storage(other.storage) {}

    T* allocate(size_t n) {
        return reinterpret_cast<T*>(storage->allocate(n * sizeof(T), alignof(T)));
    }

    void deallocate(T*, size_t) {}

    template <typename U>
    // NOLINTNEXTLINE
    struct rebind {
        using other = StackAllocator<U, N>;
    };

};
