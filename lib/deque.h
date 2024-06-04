#include <iterator>
#include <algorithm>
#include <compare>

template <typename T>
class Deque {
private:
    template <typename Value>
    // NOLINTNEXTLINE
    class base_iterator {
    public:
        using difference_type = ptrdiff_t;
        using value_type = Value;
        using reference = value_type&;
        using pointer = value_type*;
        using iterator_category = std::random_access_iterator_tag;

    private:
        // NOLINTBEGIN
        value_type** _chunk = nullptr;
        value_type* _item = nullptr;
        size_t _in_chunk_idx = 0;
        // NOLINTEND

        void _dereferenceChunk() {
            _item = *_chunk + _in_chunk_idx;
        }

    public:
        base_iterator() = default;
        base_iterator(value_type** chunk, value_type* item, size_t in_chunk_idx) :
                _chunk(chunk), _item(item), _in_chunk_idx(in_chunk_idx) {}
        base_iterator(value_type** chunk, size_t in_chunk_idx) : base_iterator(chunk, nullptr, in_chunk_idx) {
            _dereferenceChunk();
        }

        value_type& operator*() {
            return *_item;
        }

        value_type* operator->() {
            return _item;
        }

        base_iterator& operator++() {
            if (_in_chunk_idx < _chunk_size - 1) {
                ++_in_chunk_idx;
                if (_item) {
                    ++_item;
                }
                return *this;
            }
            ++_chunk;
            _in_chunk_idx = 0;
            _dereferenceChunk();
            return *this;
        }

        base_iterator operator++(int) {
            auto copy = *this;
            ++(*this);
            return copy;
        }

        base_iterator& operator--() {
            if (_in_chunk_idx > 0) {
                --_in_chunk_idx;
                if (_item) {
                    --_item;
                }
                return *this;
            }
            --_chunk;
            _in_chunk_idx = _chunk_size - 1;
            _dereferenceChunk();
            return *this;
        }

        base_iterator operator--(int) {
            auto copy = *this;
            --(*this);
            return copy;
        }

        base_iterator& operator+=(difference_type diff) {
            if (diff < 0) {
                return *this -= -diff;
            }

            auto st_diff = static_cast<size_t>(diff);
            if (_in_chunk_idx + st_diff < _chunk_size) {
                _in_chunk_idx += st_diff;
                if (_item) {
                    _item += st_diff;
                }
                return *this;
            }

            st_diff += _in_chunk_idx;
            _chunk += st_diff / _chunk_size;
            _in_chunk_idx = st_diff % _chunk_size;
            _dereferenceChunk();
            return *this;
        }

        base_iterator& operator-=(difference_type diff) {
            if (diff < 0) {
                return *this += -diff;
            }

            auto st_diff = static_cast<size_t>(diff);
            if (_in_chunk_idx >= st_diff) {
                _in_chunk_idx -= st_diff;
                if (_item) {
                    _item -= st_diff;
                }
                return *this;
            }

            st_diff += _chunk_size - _in_chunk_idx - 1;
            _chunk -= st_diff / _chunk_size;
            _in_chunk_idx = _chunk_size - st_diff % _chunk_size - 1;
            _dereferenceChunk();
            return *this;
        }

        base_iterator operator+(difference_type diff) const {
            auto copy = *this;
            copy += diff;
            return copy;
        }

        base_iterator operator-(difference_type diff) const {
            auto copy = *this;
            copy -= diff;
            return copy;
        }

        difference_type operator-(const base_iterator<value_type>& other) const {
            if (_chunk == other._chunk) {
                return static_cast<difference_type>(_in_chunk_idx - other._in_chunk_idx);
            }
            return static_cast<difference_type>(static_cast<size_t>(_chunk - other._chunk) * _chunk_size + _in_chunk_idx - other._in_chunk_idx);
        }

        std::strong_ordering operator<=>(const base_iterator<value_type>& other) const {
            auto chunks_ordering = _chunk <=> other._chunk;
            if (chunks_ordering != std::strong_ordering::equal) {
                return chunks_ordering;
            }
            return _in_chunk_idx <=> other._in_chunk_idx;
        };

        bool operator==(const base_iterator<value_type>& other) const {
            return ((*this) <=> other) == std::strong_ordering::equal;
        }

        operator base_iterator<const value_type>() const {
            return base_iterator<const value_type>(const_cast<const T**>(_chunk), const_cast<const T*>(_item), _in_chunk_idx);
        }
    };

public:
    using value_type = T;
    using iterator = base_iterator<T>;
    using const_iterator = base_iterator<const T>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

private:
    // NOLINTBEGIN
    static const size_t _chunk_size = 32;

    size_t _chunks_count = 0;
    T** _chunks_begin = nullptr;

    iterator _begin = iterator();
    iterator _end = iterator();

    // NOLINTEND

    static T** _allocate(size_t chunks_count) {
        return new T*[chunks_count];
    }

    static void _deallocate(T** chunks) {
        delete[] chunks;
    }

    static void _deallocateChunks(T** chunks, size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            delete[] reinterpret_cast<char*>(chunks[i]);
        }
    }

    static void _allocateChunks(T** chunks, size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            try {
                chunks[i] = reinterpret_cast<T*>(new char[sizeof(T) * _chunk_size]);
            } catch (std::bad_alloc&) {
                _deallocateChunks(chunks, begin, i);
                _deallocate(chunks);
                throw;
            }
        }
    }

    static void _destroyItems(iterator begin, iterator end) {
        for (auto it = begin; it != end; ++it) {
            it->~T();
        }
    }

    iterator _chunksBeginIterator() const {
        return iterator(_chunks_begin, 0);
    }

    iterator _chunksEndIterator() const {
        return iterator(_chunks_begin + _chunks_count - 1, _chunk_size - 1);
    }

    ptrdiff_t _diffFromChunksBeginIterator(const iterator& it) const {
        return it - _chunksBeginIterator();
    }

    ptrdiff_t _diffBetweenChunksBeginAndBegin() const {
        return _diffFromChunksBeginIterator(_begin);
    }

    ptrdiff_t _diffBetweenChunksBeginAndEnd() const {
        return _diffFromChunksBeginIterator(_end);
    }

    void _swap(Deque& other) noexcept {
        std::swap(_chunks_begin, other._chunks_begin);
        std::swap(_chunks_count, other._chunks_count);
        std::swap(_begin, other._begin);
        std::swap(_end, other._end);
    }

    size_t _getChunkIdx(const iterator& it) const {
        return _diffFromChunksBeginIterator(it) / _chunk_size;
    }

    std::tuple<T**, size_t, size_t> _increaseChunksCountAllocateAndCopy() {
        size_t unused_chunks_left = static_cast<size_t>(_diffFromChunksBeginIterator(_begin)) / _chunk_size;
        size_t unused_chunks_right = _chunks_count - static_cast<size_t>(_diffFromChunksBeginIterator(_end) - 1) / _chunk_size - 1;

        size_t add_chunks_count = _chunks_count;
        if ((add_chunks_count + unused_chunks_right - unused_chunks_left) % 2) {
            ++add_chunks_count;
        }

        size_t add_chunks_left = (add_chunks_count + unused_chunks_right - unused_chunks_left) / 2;

        size_t new_chunks_count = _chunks_count + add_chunks_count;

        T** new_chunks = _allocate(new_chunks_count);
        _allocateChunks(new_chunks, 0, add_chunks_left);
        try {
            _allocateChunks(new_chunks, add_chunks_left + _chunks_count, new_chunks_count);
        } catch (...) {
            _deallocateChunks(new_chunks, 0, add_chunks_left);
            _deallocate(new_chunks);
            throw;
        }

        try {
            std::copy(_chunks_begin, _chunks_begin + _chunks_count, new_chunks + add_chunks_left);
        } catch (...) {
            _deallocateIncreaseChunks(new_chunks, new_chunks_count, add_chunks_left);
            throw;
        }


        return {new_chunks, new_chunks_count, add_chunks_left};
    }

    void _deallocateIncreaseChunks(T** new_chunks, size_t new_chunks_count, size_t add_chunks_left) {
        _deallocateChunks(new_chunks, 0, add_chunks_left);
        _deallocateChunks(new_chunks, add_chunks_left + _chunks_count, new_chunks_count);
        _deallocate(new_chunks);
    }

    void _increaseChunksCountDeallocate(T** new_chunks, size_t new_chunks_count, size_t add_chunks_left) {
        ptrdiff_t begin_shift = _diffBetweenChunksBeginAndBegin() + static_cast<ptrdiff_t>(_chunk_size * add_chunks_left);
        ptrdiff_t end_shift = _diffBetweenChunksBeginAndEnd() + static_cast<ptrdiff_t>(_chunk_size * add_chunks_left);

        _deallocate(_chunks_begin);

        _chunks_begin = new_chunks;
        _chunks_count = new_chunks_count;

        _begin = _chunksBeginIterator() + begin_shift;
        _end = _chunksBeginIterator() + end_shift;
    }

    void _increaseAndInsert(iterator pos, const T& item) {
        ptrdiff_t from_old_begin_to_pos_shift = pos - _begin;
        auto [new_chunks, new_chunks_count, add_chunks_left] = _increaseChunksCountAllocateAndCopy();

        ptrdiff_t begin_shift = _diffBetweenChunksBeginAndBegin() + static_cast<ptrdiff_t>(_chunk_size * add_chunks_left);
        pos = iterator(new_chunks, 0) + begin_shift + from_old_begin_to_pos_shift;
        try {
            new (&(*pos)) T(item);
        } catch (...) {
            _deallocateIncreaseChunks(new_chunks, new_chunks_count, add_chunks_left);
            throw;
        }
        _increaseChunksCountDeallocate(new_chunks, new_chunks_count, add_chunks_left);
    }

public:
    Deque() : _chunks_count(1), _chunks_begin(_allocate(_chunks_count))
    {
        _allocateChunks(_chunks_begin, 0, _chunks_count);
        _begin = _end = _chunksBeginIterator() + _chunk_size / 2;
    }

    explicit Deque(size_t size)
            : _chunks_count(((size + 1) + 1) / _chunk_size + 1),
              _chunks_begin(_allocate(_chunks_count)),
              _begin(iterator()), _end(iterator()) {
        _construct_several_items(size);
    }

    Deque(size_t size, const T& item)
            : _chunks_count(((size + 1) + 1) / _chunk_size + 1),
              _chunks_begin(_allocate(_chunks_count)),
              _begin(iterator()), _end(iterator()) {
        _construct_several_items(size, item);
    }

    template <typename ...Args>
    void _construct_several_items(size_t size, const Args&... args) {
        _allocateChunks(_chunks_begin, 0, _chunks_count);
        _begin = _chunksBeginIterator() + 1;
        _end = _begin + static_cast<ptrdiff_t>(size);

        for (auto it = _begin; it != _end; ++it) {
            try {
                new (&(*it)) T(args...);
            } catch (...) {
                _destroyItems(_begin, it);
                _deallocateChunks(_chunks_begin, 0, _chunks_count);
                _deallocate(_chunks_begin);
                throw;
            }
        }
    }

    Deque(const Deque& other) : _chunks_count(other._chunks_count), _chunks_begin(_allocate(_chunks_count)) {
        _allocateChunks(_chunks_begin, 0, _chunks_count);

        ptrdiff_t begin_shift = other._diffBetweenChunksBeginAndBegin();
        ptrdiff_t end_shift = other._diffBetweenChunksBeginAndEnd();

        _begin = _chunksBeginIterator() + begin_shift;
        _end = _chunksBeginIterator() + end_shift;

        for (auto it = _begin, it_other = other._begin; it != _end; ++it, ++it_other) {
            try {
                new (&(*it)) T(*it_other);
            } catch (...) {
                _destroyItems(_begin, it);
                _deallocateChunks(_chunks_begin, 0, _chunks_count);
                _deallocate(_chunks_begin);
                throw;
            }
        }
    }

    Deque& operator=(const Deque& other) {
        if (this != &other) {
            Deque copy(other);
            _swap(copy);
        }
        return *this;
    }

    ~Deque() {
        for (auto it = _begin; it != _end; ++it) {
            it->~T();
        }
        for (size_t i = 0; i < _chunks_count; ++i) {
            delete[] reinterpret_cast<char*>(_chunks_begin[i]);
        }
        delete[] _chunks_begin;
    }

    size_t size() const {
        return static_cast<size_t>(_end - _begin);
    }

    T& operator[](size_t index) {
        return *(_begin + static_cast<ptrdiff_t>(index));
    }

    const T& operator[](size_t index) const {
        return *(_begin + static_cast<ptrdiff_t>(index));
    }

    T& at(size_t index) {
        if (index >= size()) {
            throw std::out_of_range("Index is out of range");
        }
        return this->operator[](index);
    }

    const T& at(size_t index) const {
        if (index >= size()) {
            throw std::out_of_range("Index is out of range");
        }
        return this->operator[](index);
    }

    iterator begin() {
        return _begin;
    }

    const_iterator begin() const {
        return _begin;
    }

    iterator end() {
        return _end;
    }

    const_iterator end() const {
        return _end;
    }

    const_iterator cbegin() const {
        return _begin;
    }

    const_iterator cend() const {
        return _end;
    }

    reverse_iterator rbegin() {
        return std::reverse_iterator(end());
    }

    const_reverse_iterator rbegin() const {
        return std::reverse_iterator(end());
    }

    reverse_iterator rend() {
        return std::reverse_iterator(begin());
    }

    const_reverse_iterator rend() const {
        return std::reverse_iterator(begin());
    }

    const_reverse_iterator crbegin() const {
        return rbegin();
    }

    const_reverse_iterator crend() const {
        return rend();
    }


    void push_back(const T& item) {
        if (_end + 1 == _chunksEndIterator()) {
            _increaseAndInsert(_end, item);
        } else {
            new (&(*_end)) T(item);
        }
        ++_end;
    }

    void push_front(const T& item) {
        if (_begin - 1 == _chunksBeginIterator()) {
            _increaseAndInsert(_begin - 1, item);
        } else {
            new (&(*(_begin - 1))) T(item);
        }
        --_begin;
    }

    void pop_back() {
        (_end - 1)->~T();
        --_end;
    }

    void pop_front() {
        _begin->~T();
        ++_begin;
    }

    iterator insert(iterator insert_it, const T& item) {
        ptrdiff_t shift = insert_it - _begin;
        push_back(item);
        insert_it = _begin + shift;
        for (auto it = _end - 1; it != insert_it; --it) {
            std::swap(*it, *(it - 1));
        }
        return insert_it;
    }

    iterator erase(iterator erase_it) {
        erase_it->~T();
        for (auto it = erase_it + 1; it != _end; ++it) {
            new (&(*(it - 1))) T(*it);
            it->~T();
        }
        --_end;
        return erase_it;
    }
};