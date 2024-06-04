#include <functional>
#include <cmath>
#include <vector>
#include <stdexcept>

template<typename Key, typename Value, typename Hash = std::hash<Key>,
        typename EqualTo = std::equal_to<Key>,
        typename Alloc = std::allocator<std::pair<const Key, Value>>>
class UnorderedMap {
public:
    using NodeType = std::pair<const Key, Value>;
private:
    template<typename T, typename ListAlloc = std::allocator<T> >
    class List {
        friend UnorderedMap<Key, Value, Hash, EqualTo, Alloc>;
    private:
        struct BaseNode {
            BaseNode* prev = nullptr;
            BaseNode* next = nullptr;
        };

        struct Node : BaseNode {
            T item;
            size_t hash_t = 0;
        };

        template<typename IteratorT>
        class BaseIterator {
        public:
            using difference_type = ptrdiff_t;
            using value_type = IteratorT;
            using reference = value_type&;
            using pointer = value_type*;
            using iterator_category = std::bidirectional_iterator_tag;

        private:
            BaseNode* node = nullptr;

        public:
            BaseIterator() = default;

            explicit BaseIterator(BaseNode* node) : node(node) {}

            value_type& operator*() const {
                return static_cast<Node*>(node)->item;
            }

            value_type* operator->() const {
                return &operator*();
            }

            BaseIterator& operator++() {
                node = node->next;
                return *this;
            }

            BaseIterator operator++(int) {
                auto copy = *this;
                ++(*this);
                return copy;
            }

            BaseIterator& operator--() {
                node = node->prev;
                return *this;
            }

            BaseIterator operator--(int) {
                auto copy = *this;
                --(*this);
                return copy;
            }


            bool operator==(const BaseIterator<value_type>& other) const = default;

            operator BaseIterator<const value_type>() const {
                return BaseIterator<const value_type>(node);
            }

            BaseNode* getNode() {
                return node;
            }

            bool is_null() const {
                return !node;
            }
        };

    public:
        using value_type = T;
        using iterator = BaseIterator<T>;
        using const_iterator = BaseIterator<const T>;

    private:
        using NodeAlloc = typename std::allocator_traits<ListAlloc>::template rebind_alloc<Node>;
        using NodeTraits = std::allocator_traits<NodeAlloc>;
        [[no_unique_address]] NodeAlloc node_alloc;
        BaseNode fake;
        size_t container_size = 0;

        void _install_node_before(BaseNode* place, Node* node) {
            node->next = place;
            node->prev = place->prev;

            place->prev = place->prev->next = node;

            ++container_size;
        }

        void _install_back_node(Node* node) {
            _install_node_before(&fake, node);
        }

        void _install_front_node(Node* node) {
            _install_node_before(fake.next, node);
        }

        BaseNode* _extract_node(Node* node) {
            node->prev->next = node->next;
            node->next->prev = node->prev;

            BaseNode* next_node = node->next;
            node->next = node->prev = nullptr;
            --container_size;
            return next_node;
        }

        Node* _extract_back_node() {
            Node* node = static_cast<Node*>(fake.prev);
            _extract_node(node);
            return node;
        }

        void _destroy_node(Node* node) {
            NodeTraits::destroy(node_alloc, node);
            NodeTraits::deallocate(node_alloc, node, 1);
        }

        void _clear() {
            while (container_size) {
                _destroy_node(_extract_back_node());
            }
        }

        template<typename... Args>
        Node* _create_node(Args&& ... args) {
            Node* node = NodeTraits::allocate(node_alloc, 1);
            try {
                NodeTraits::construct(node_alloc, &node->Node::item, std::forward<Args>(args)...);
                node->next = node->prev = nullptr;
                node->hash_t = 0;
            } catch (...) {
                NodeTraits::deallocate(node_alloc, node, 1);
                throw;
            }
            return node;
        }

        template<typename... Args>
        Node* _create_node_with_hash(size_t hash_t, Args&&... args) {
            Node* node = _create_node(std::forward<Args>(args)...);
            node->hash_t = hash_t;
            return node;
        }

        void _after_copy_constructed(const List& other) {
            try {
                for (auto it = other.begin(); it != other.end(); ++it) {
                    auto node = static_cast<Node*>(it.getNode());
                    _install_back_node(_create_node_with_hash(node->hash_t, node->item));
                }
            } catch (...) {
                _clear();
                throw;
            }
        }



    public:
        List() : List(ListAlloc()) {}

        explicit List(const ListAlloc& alloc) : node_alloc(alloc), fake(), container_size(0) {
            fake.next = fake.prev = &fake;
        };

        List(const List& other) : List(std::allocator_traits<NodeAlloc>::select_on_container_copy_construction(
                other.node_alloc)) {
            _after_copy_constructed(other);
        }

        List(const List& other, const ListAlloc& alloc) : List(alloc) {
            _after_copy_constructed(other);
        }

        List(List&& other) noexcept : node_alloc(std::move(other.node_alloc)) {
            _after_move_constructed(std::move(other));
        }
        List(List&& other, const ListAlloc& alloc) noexcept : node_alloc(alloc) {
            _after_move_constructed(std::move(other));
        }

        void _after_move_constructed(List&& other) noexcept {
            if (other.container_size) {
                fake.next = other.fake.next;
                fake.prev = other.fake.prev;
                fake.next->prev = fake.prev->next = &fake;
                other.fake.next = other.fake.prev = &other.fake;
            } else {
                fake.next = fake.prev = &fake;
            }


            container_size = other.container_size;
            other.container_size = 0;
        }

        List& operator=(const List& other) {
            if (this != &other) {
                List copy(other, std::allocator_traits<ListAlloc>::propagate_on_container_copy_assignment::value
                                 ? other.node_alloc : node_alloc);
                swap(copy);
            }
            return *this;
        }

        void _move_assign(List&& other, std::true_type) {
            if (std::allocator_traits<ListAlloc>::propagate_on_container_move_assignment::value) {
                node_alloc = std::move(other.node_alloc);
            }
            _after_move_constructed(std::move(other));
        }

        void _move_assign(List&& other, std::false_type) {
            if (node_alloc == other.node_alloc) {
                _move_assign(std::move(other), std::true_type());
                return;
            }
            while (other.size()) {
                Node* node = other._extract_back_node();
                push_front(node->hash_t, std::move(node->item));
                other._destroy_node(node);
            }
        }

        List& operator=(List&& other) noexcept {
            _clear();
            _move_assign(std::move(other), typename std::allocator_traits<ListAlloc>::propagate_on_container_move_assignment());
            return *this;
        }

        void swap(List& other) {
            if (container_size) {
                fake.next->prev = fake.prev->next = &other.fake;
            } else {
                fake.next = fake.prev = &other.fake;
            }

            if (other.container_size) {
                other.fake.next->prev = other.fake.prev->next = &fake;
            } else {
                other.fake.next = other.fake.prev = &fake;
            }

            std::swap(fake, other.fake);
            std::swap(container_size, other.container_size);

            if (std::allocator_traits<ListAlloc>::propagate_on_container_swap::value) {
                std::swap(node_alloc, other.node_alloc);
            }
        }

        NodeAlloc get_allocator() const {
            return node_alloc;
        }

        size_t size() const {
            return container_size;
        }

        iterator begin() {
            return iterator(fake.next);
        }

        const_iterator begin() const {
            return iterator(fake.next);
        }

        iterator end() {
            return iterator(&fake);
        }

        const_iterator end() const {
            return iterator(const_cast<List::BaseNode*>(&fake));
        }

        const_iterator cbegin() const {
            return begin();
        }

        const_iterator cend() const {
            return end();
        }

        template <typename... Args>
        iterator push_front(size_t hash_t, Args&&... args) {
            Node* node = _create_node_with_hash(hash_t, std::forward<Args>(args)...);
            return push_front(hash_t, node);
        }

        iterator push_front(size_t, Node* node) {
            _install_front_node(node);
            return iterator(node);
        }

        template <typename... Args>
        iterator insert(const_iterator it, size_t hash_t, Args&&... args) {
            Node* node = _create_node_with_hash(hash_t, std::forward<Args>(args)...);
            return insert(it, hash_t, node);
        }

        iterator insert(const_iterator it, size_t, Node* node) {
            _install_node_before(it.getNode(), node);
            return iterator(node);
        }

        iterator erase(const_iterator it) {
            return iterator(_extract_node(static_cast<Node*>(it.getNode())));
        }

        ~List() {
            _clear();
        }
    };

    // NOLINTBEGIN
    static const size_t INIT_BUCKET_SIZE = 2;
    constexpr static const float INIT_MAX_LOAD_FACTOR = 1.0;
    // NOLINTEND

    using ListT = List<NodeType, Alloc>;
    using ListBaseNode =  typename ListT::BaseNode;
    using ListNode = typename ListT::Node;

    using BucketsVectorAllocT = typename std::allocator_traits<Alloc>::template rebind_alloc<typename ListT::iterator>;
    using ListNodePVectorAllocT = typename std::allocator_traits<Alloc>::template rebind_alloc<ListNode*>;

    // NOLINTBEGIN
    [[ no_unique_address ]] Hash _hash;
    [[ no_unique_address ]] EqualTo _equal;
    [[ no_unique_address ]] Alloc _alloc;
    ListT _items;
    std::vector<typename ListT::iterator, BucketsVectorAllocT> _buckets;
    float _max_load_factor = INIT_MAX_LOAD_FACTOR;
    // NOLINTEND

public:
    using iterator = typename ListT::iterator;
    using const_iterator = typename ListT::const_iterator;

    UnorderedMap() : UnorderedMap(Alloc()) {}

    explicit UnorderedMap(const Alloc& alloc) : UnorderedMap(INIT_BUCKET_SIZE, alloc) {}

    explicit UnorderedMap(size_t bucket_count,
                          const Hash& hash = Hash(),
                          const EqualTo& equal = EqualTo(),
                          const Alloc& alloc = Alloc())
            : _hash(hash), _equal(equal), _alloc(alloc),
              _items(_alloc),
              _buckets(std::max(bucket_count, static_cast<size_t>(1)), _alloc) {
    }

    UnorderedMap(size_t bucket_count, const Alloc& alloc)
            : UnorderedMap(bucket_count, Hash(), EqualTo(), alloc) {}

    UnorderedMap(size_t bucket_count,
                 const Hash& hash,
                 const Alloc& alloc)
            : UnorderedMap(bucket_count, hash, EqualTo(), alloc) {}

    UnorderedMap(const UnorderedMap& other) :
            UnorderedMap(other, std::allocator_traits<Alloc>::select_on_container_copy_construction(other._alloc)) {}


    UnorderedMap(const UnorderedMap& other, const Alloc& alloc)
            : _hash(other._hash), _equal(other._equal), _alloc(alloc),
              _items(other._items, alloc), _buckets(other._buckets.size(), alloc),
              _max_load_factor(other._max_load_factor) {
        _restore_buckets();
    }

    UnorderedMap(UnorderedMap&& other) noexcept
            : _hash(std::move(other._hash)), _equal(std::move(other._equal)),
              _alloc(std::move(other._alloc)), _items(std::move(other._items)), _buckets(std::move(other._buckets)),
              _max_load_factor(other._max_load_factor) {
        _after_move_constructed(std::move(other));
    }

    UnorderedMap(UnorderedMap&& other, const Alloc& alloc) noexcept
            : _hash(std::move(other._hash)), _equal(std::move(other._equal)),
              _alloc(alloc), _items(std::move(other._items), alloc), _buckets(std::move(other._buckets), alloc),
              _max_load_factor(other._max_load_factor) {
        _after_move_constructed(std::move(other));
    }

    void swap(UnorderedMap& other) {
        std::swap(_hash, other._hash);
        std::swap(_equal, other._equal);
        std::swap(_buckets, other._buckets);
        std::swap(_max_load_factor, other._max_load_factor);

        _items.swap(other._items);
        if (std::allocator_traits<Alloc>::propagate_on_container_swap::value) {
            std::swap(_alloc, other._alloc);
        }
    }

    UnorderedMap& operator=(const UnorderedMap& other) {
        if (this != &other) {
            UnorderedMap copy(other, std::allocator_traits<Alloc>::propagate_on_container_copy_assignment::value
                                     ? other._alloc : _alloc);
            swap(copy);
        }
        return *this;
    }

    UnorderedMap& operator=(UnorderedMap&& other) noexcept {
        _hash = std::move(other._hash);
        _equal = std::move(other._equal);
        _items = std::move(other._items);
        _buckets = std::move(other._buckets);
        _max_load_factor = other._max_load_factor;
        _after_move_constructed(std::move(other));

        if (!std::allocator_traits<Alloc>::propagate_on_container_move_assignment::value && (_alloc != other._alloc)) {
            _restore_buckets();
        } else if (std::allocator_traits<Alloc>::propagate_on_container_move_assignment::value) {
            _alloc = std::move(other._alloc);
        }
        return *this;
    }

    ~UnorderedMap() = default;

    void rehash(size_t new_bucket_count) {
        new_bucket_count = std::max(new_bucket_count, _get_buckets_for_max_load_factor(size()));
        if (new_bucket_count == _buckets.size()) {
            return;
        }
        std::vector<ListNode*, ListNodePVectorAllocT> nodes(_alloc);
        nodes.reserve(size());
        while (_items.container_size) {
            nodes.push_back(_items._extract_back_node());
        }

        _buckets.assign(new_bucket_count, iterator());

        for (auto node : nodes) {
            _insert_node(node);
        }
    }

    Value& operator[](const Key& key) {
        return _iterator_get_node(_insert(key,
                                          std::piecewise_construct,
                                          std::forward_as_tuple(key),
                                          std::forward_as_tuple()).first)->item.second;
    }

    Value& operator[](Key&& key) {
        return _iterator_get_node(_insert(key,
                                          std::piecewise_construct,
                                          std::forward_as_tuple(std::forward<Key>(key)),
                                          std::forward_as_tuple()).first)->item.second;
    }

    Value& at(const Key& key) {
        return const_cast<Value&>(_at(key));
    }

    const Value& at(const Key& key) const {
        return _at(key);
    }

    std::pair<iterator, bool> insert(const NodeType& value) {
        return _insert(value.first, value);
    }

    std::pair<iterator, bool> insert(NodeType&& value) {
        return _insert(value.first, std::forward<NodeType>(value));
    }

    template <typename Pair>
    std::__enable_if_t<std::is_constructible<NodeType, Pair&&>::value,
            std::pair<iterator, bool>>
    insert(Pair&& value) {
        return _insert(value.first, std::forward<Pair>(value));
    }

    template <class InputIt>
    void insert(InputIt first, InputIt last) {
        std::vector<ListNode*, ListNodePVectorAllocT> nodes(_alloc);
        try {
            while (first != last) {
                nodes.push_back(_create_node_with_hash(*first));
                ++first;
            }
        } catch (...) {
            for (auto node : nodes) {
                _items._destroy_node(node);
            }
            throw;
        }

        for (auto node : nodes) {
            _insert_node(node);
        }
    }

    template <typename ...Args>
    std::pair<iterator, bool> emplace(Args&&... args) {
        return _insert_node(_create_node_with_hash(std::forward<Args>(args)...));
    }

    iterator erase(const_iterator const_pos) {
        iterator pos = _iterator_const_cast(const_pos);
        size_t bucket = _get_node_bucket(_iterator_get_node(pos));
        if (_buckets[bucket] == pos) {
            auto next_it = pos;
            ++next_it;
            if ((next_it == _items.end()) || (_get_node_bucket(_iterator_get_node(next_it)) != bucket)) {
                _buckets[bucket] = iterator();
            } else {
                _buckets[bucket] = next_it;
            }
        }

        return _items.erase(pos);
    }

    iterator erase(const_iterator first, const_iterator last) {
        while (first != last) {
            first = erase(first);
        }
        return _iterator_const_cast(first);
    }

    iterator find(const Key& key) {
        return _iterator_const_cast(_find(key));
    }

    const_iterator find(const Key& key) const {
        return _find(key);
    }

    void reserve(size_t count) {
        rehash(_get_buckets_for_max_load_factor(count));
    }

    Alloc get_allocator() const noexcept {
        return _items.get_allocator();
    }

    size_t bucket_count() const {
        return _buckets.size();
    }

    size_t size() const {
        return _items.size();
    }

    float max_load_factor() const {
        return _max_load_factor;
    }

    void max_load_factor(float ml) {
        _max_load_factor = ml;
        if (_is_rehash_required()) {
            rehash(_get_buckets_for_max_load_factor(size()));
        }
    }

    float load_factor() const {
        return static_cast<float>(size()) / static_cast<float>(_buckets.size());
    }

    iterator begin() {
        return _items.begin();
    }

    const_iterator begin() const {
        return _items.begin();
    }

    iterator end() {
        return _items.end();
    }

    const_iterator end() const {
        return _items.end();
    }

    const_iterator cbegin() const {
        return begin();
    }

    const_iterator cend() const {
        return end();
    }

private:
    void _restore_buckets() {
        size_t previous_bucket = _buckets.size();
        for (auto it = _items.begin(); it != _items.end(); ++it) {
            auto current_bucket = _iterator_get_node(it)->hash_t % _buckets.size();
            if (current_bucket != previous_bucket) {
                _buckets[current_bucket] = it;
                previous_bucket = current_bucket;
            }
        }
    }

    void _after_move_constructed(UnorderedMap&& other) {
        other._buckets.assign(1, iterator());
        other._max_load_factor = INIT_MAX_LOAD_FACTOR;
    }

    void _check_rehash() {
        if (_is_rehash_required()) {
            // NOLINTNEXTLINE
            rehash(_buckets.size() * 3 + 17);
        }
    }

    const Value& _at(const Key& key) const {
        auto [it, is_not_exist] = _find_place(key, _hash(key));
        if (!is_not_exist) {
            return _iterator_get_node(it)->item.second;
        }
        throw std::out_of_range("No such key");
    }

    template <typename KeyRef, typename ...Args>
    std::pair<iterator, bool> _insert(KeyRef&& key, Args&&... args) {
        return _insert_with_hash(_hash(key), std::forward<KeyRef>(key), std::forward<Args>(args)...);
    }

    template <typename KeyRef, typename ...Args>
    std::pair<iterator, bool> _insert_with_hash(size_t hash_t, KeyRef&& key, Args&&... args) {
        auto [const_insert_it, is_not_exist] = _find_place(key, hash_t);
        iterator insert_it = _iterator_const_cast(const_insert_it);
        if (!is_not_exist) {
            return {insert_it, false};
        }

        iterator res_it;
        if (insert_it == _items.begin()) {
            res_it = _items.push_front(hash_t, std::forward<Args>(args)...);
            _buckets[_get_hash_bucket(hash_t)] = res_it;
        } else {
            res_it = _items.insert(insert_it, hash_t, std::forward<Args>(args)...);
        }

        _check_rehash();
        return {res_it, true};
    }

    std::pair<iterator, bool> _insert_node(ListNode* node) {
        auto res = _insert_with_hash(node->hash_t, node->item.first, node);;
        if (!res.second) {
            _items._destroy_node(node);
        }
        return res;
    }

    size_t _get_hash_bucket(size_t hash_t) const {
        return hash_t % _buckets.size();
    }

    size_t _get_buckets_for_max_load_factor(size_t count) const {
        return static_cast<size_t>(std::ceil(static_cast<float>(count) / _max_load_factor));
    }

    bool _is_rehash_required() const {
        return load_factor() > _max_load_factor;
    }

    size_t _get_node_bucket(ListNode* node) const {
        return _get_hash_bucket(node->hash_t);
    }

    template <typename ...Args>
    ListNode* _create_node_with_hash(Args&&... args) {
        ListNode* node = _items._create_node(std::forward<Args>(args)...);
        node->hash_t = _hash(node->item.first);
        return node;
    }

    const_iterator _find(const Key& key) const {
        size_t hash_t = _hash(key);
        auto [it, is_not_exist] = _find_place(key, hash_t);
        if (is_not_exist) {
            return _items.cend();
        }
        return it;
    }

    std::pair<const_iterator, bool> _find_place(const Key& key, size_t hash_t) const {
        size_t bucket = _get_hash_bucket(hash_t);
        auto bucket_it = _buckets[bucket];
        if (bucket_it.is_null()) {
            return {_items.begin(), true};
        }

        while (true) {
            auto node = _iterator_get_node(bucket_it);
            if (_equal(key, node->item.first)) {
                return {bucket_it, false};
            }
            ++bucket_it;
            if ((bucket_it == _items.end()) || (_get_node_bucket(_iterator_get_node(bucket_it)) != bucket)) {
                break;
            }
        }
        return {bucket_it, true};
    }

    static iterator _iterator_const_cast(const_iterator it) {
        return iterator(const_cast<ListBaseNode*>(it.getNode()));
    }

    static ListNode* _iterator_get_node(iterator it) {
        return static_cast<ListNode*>(it.getNode());
    }

    static ListNode* _iterator_get_node(const_iterator it) {
        return static_cast<ListNode*>(const_cast<ListBaseNode*>(it.getNode()));
    }
};