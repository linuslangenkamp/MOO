#ifndef OPT_OBSERVER_H
#define OPT_OBSERVER_H

template<typename T>
class Observer {
public:
    Observer() = default;

    Observer(std::nullptr_t) noexcept : ptr(nullptr) {}

    explicit Observer(const T* p) noexcept : ptr(p) {}

    Observer& operator=(const T* p) noexcept {
        ptr = p;
        return *this;
    }

    const T& operator*() const {
        assert(ptr && "Dereferencing null Observer!");
        return *ptr;
    }

    const T* operator->() const {
        assert(ptr && "Dereferencing null Observer!");
        return ptr;
    }

    const T* get() const noexcept { return ptr; }

    explicit operator bool() const noexcept { return ptr != nullptr; }

    bool operator==(const Observer& other) const noexcept { return ptr == other.ptr; }

    bool operator!=(const Observer& other) const noexcept { return ptr != other.ptr; }

private:
    const T* ptr = nullptr;
};

#endif // OPT_OBSERVER_H
