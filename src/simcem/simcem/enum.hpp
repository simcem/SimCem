#pragma once
#include <stator/xml.hpp>

namespace simcem {
  using namespace stator::xml;

  template<class T>
  struct Enum {
    Enum(size_t v) : _value(v) {
      if (_value > size())
	stator_throw() << "Out of range enum " << _value;
    }
      
    Enum(const std::string str) {
      for (size_t i(0); i < size(); ++i)
	if (str == T::strings[i]) {
	  _value = i;
	  return;
	}
      stator_throw() << "Unrecognised enumerated value " << str;
    }

    Enum(const Attribute xml) {
      *this = Enum(std::string(xml));
    }

    operator int() const { return _value; }
    operator std::string() const { return T::strings[_value]; }

    bool operator==(const Enum& o) const {
      return _value == o._value;
    }

    bool operator==(int o) const {
      return int(_value) == o;
    }

    bool operator==(size_t o) const {
      return _value == int(o);
    }
    
    inline static size_t size() {
      return sizeof(T::strings) / sizeof(T::strings[0]);
    }
    
    size_t _value;
  };

  template<class T>
  std::ostream& operator<<(std::ostream& os, const Enum<T>& en) {
    os << std::string(en);
    return os;
  }
}
