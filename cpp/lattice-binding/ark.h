//
//  ark.h
//  lattice-binding
//
//  Created by Timothy Palpant on 3/16/13.
//  Copyright (c) 2013 Timothy Palpant. All rights reserved.
//

#ifndef lattice_binding_ark_h
#define lattice_binding_ark_h

#include <string>

namespace config {
  class Ark {
  private:
    bool open_;
    
  public:
    static Ark load(std::string& path);
    static Ark fromArgv(const char* argv[]);
    bool has(std::string& key);
    
    template <typename T>
    T get(std::string& key);
    
    template <typename T>
    void set(std::string& key, T value);
    
    void update(Ark& ark);
    
    size_t size();
    
    bool is_open() const { return open_; }
    
    void set_open(bool open) { open_ = open; }
    
    std::string to_string();
  };
}

#endif
