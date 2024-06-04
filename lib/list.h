#include "stackallocator.h"

template <typename T, typename Alloc = std::allocator<T> >
class List {
private:
    struct BaseNode {
        BaseNode* prev = nullptr;
        BaseNode* next = nullptr;
    };

    struct Node : BaseNode {
        T item;

        Node() = default;
        explicit Node(const T& item) : BaseNode(), item(item) {}
    };

    template <typename Value>
    class BaseIterator {
    public:
        using difference_type = ptrdiff_t;
        using value_type = Value;
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
    };

public:
    using value_type = T;
    using iterator = BaseIterator<T>;
    using const_iterator = BaseIterator<const T>;
    using reverse_iterator = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

private:
    using NodeAlloc = typename std::allocator_traits<Alloc>::template rebind_alloc<Node>;
    using NodeTraits = std::allocator_traits<NodeAlloc>;
    [[ no_unique_address ]] NodeAlloc node_alloc;
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

    Node* _extract_front_node() {
        Node* node = static_cast<Node *>(fake.next);
        _extract_node(node);
        return node;
    }

    Node* _extract_back_node() {
        Node* node = static_cast<Node *>(fake.prev);
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

    template <typename... Args>
    Node* _create_node(const Args&... args) {
        Node* node = NodeTraits::allocate(node_alloc, 1);
        try {
            NodeTraits::construct(node_alloc, node, args...);
        } catch (...) {
            NodeTraits::deallocate(node_alloc, node, 1);
            throw;
        }
        return node;
    }

    template <typename... Args>
    void _constructListWithSeveralObjects(size_t n, const Args&... args) {
        try {
            for (size_t i = 0; i < n; ++i) {
                Node* node = _create_node(args...);
                _install_back_node(node);
            }
        } catch (...) {
            _clear();
            throw;
        }
    }

    void _swap(List& other) {
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
        std::swap(node_alloc, other.node_alloc);
        std::swap(container_size, other.container_size);
    }

public:
    List() : List(Alloc()) {}

    explicit List(size_t n) : List(n, Alloc()) {}

    List(size_t n, const T& item) : List(n, item, Alloc()) {}

    explicit List(const Alloc& alloc) : node_alloc(alloc), fake(), container_size(0) {
        fake.next = fake.prev = &fake;
    };

    List(size_t n, const Alloc& alloc) : List(alloc) {
        _constructListWithSeveralObjects(n);
    }

    List(size_t n, const T& item, const Alloc& alloc) : List(alloc) {
        _constructListWithSeveralObjects(n, item);
    }

    List(const List& other, const Alloc& alloc) : List(alloc) {
        try {
            for (const T& item : other) {
                _install_back_node(_create_node(item));
            }
        } catch (...) {
            _clear();
            throw;
        }
    }

    List(const List& other) : List(other, std::allocator_traits<Alloc>::select_on_container_copy_construction(other.node_alloc)) {}

    List& operator=(const List& other) {
        if (this != &other) {
            List copy(other, std::allocator_traits<Alloc>::propagate_on_container_copy_assignment::value ? other.node_alloc : node_alloc);
            _swap(copy);
        }
        return *this;
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
        _install_back_node(_create_node(item));
    }

    void push_front(const T& item) {
        _install_front_node(_create_node(item));
    }

    void pop_back() {
        _destroy_node(_extract_back_node());
    }

    void pop_front() {
        _destroy_node(_extract_front_node());
    }

    iterator insert(const_iterator it, const T& item) {
        Node* node = _create_node(item);
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