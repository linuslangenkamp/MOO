#ifndef OPT_OBSERVER_PTR_H
#define OPT_OBSERVER_PTR_H

#include <vector>
#include <memory>

template<typename T> 
class observable;

template<typename T>
class observer_ptr {
public:
    observer_ptr() = default;

    observer_ptr(std::nullptr_t) noexcept : ptr(nullptr) {}

    const T& operator*() const {
        assert(ptr && "Dereferencing null observer_ptr!");
        return **ptr;
    }

    const T* operator->() const {
        assert(ptr && "Dereferencing null observer_ptr!");
        return *ptr;
    }

    const T* get() const noexcept { return *ptr; }

    explicit operator bool() const noexcept { return ptr.get() != nullptr && *ptr != nullptr; }

    bool operator==(const observer_ptr& other) const noexcept { return ptr == other.ptr; }

    bool operator!=(const observer_ptr& other) const noexcept { return ptr != other.ptr; }

private:
    friend class observable<T>;

    std::shared_ptr<const T*> ptr = nullptr;
};

template<typename T> 
class observable {
public:
    observable() {
        observed_this = std::make_shared<const T*>(reinterpret_cast<const T*>(this));
    }

    virtual ~observable() {
        *observed_this = nullptr;
    }

    observer_ptr<T> get_observer() const {
         observer_ptr<T> observer;
         observer.ptr = observed_this;

        return observer;
    }

private: 
    std::shared_ptr<const T*> observed_this;
};


#endif // OPT_OBSERVER_PTR_H
