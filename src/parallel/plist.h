//
// BAGEL - Parallel electron correlation program.
// Filename: plist.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#ifndef __SRC_PARALLEL_PLIST_H
#define __SRC_PARALLEL_PLIST_H

#include <stddef.h>
#include <memory>
#include <utility>
#include <iostream>
#include <thread>
#include <vector>

#include <src/util/constants.h>
#include <src/parallel/mpi_interface.h>

#include <bagel_config.h>
#include <src/parallel/resources.h>

namespace bagel {

// A very simple threadsafe singly-linked list modeled after the one described in
// 'C++ concurrency in action' by Anthony Williams

template <typename T>
class ParallelList {
  protected:
    struct Node {
      mutable std::mutex m;
      std::shared_ptr<T> data;

      std::unique_ptr<Node> next;

      Node() : next() {}
      Node(std::shared_ptr<T> v) : data(v), next() {}
    };

    Node head_;

  public:
    ParallelList() { }
    virtual ~ParallelList() { remove_if([] (const T& t) { return true;}); }

    ParallelList(const ParallelList& o)=delete;
    ParallelList& operator=(const ParallelList& o)=delete;

    bool empty() const {
      std::lock_guard<std::mutex> lock(head_.m);
      return (!head_.next);
    }    

    template <typename... Args>
    void emplace_front(Args&&... args) {
      push_front(std::make_shared<T>(std::forward<Args>(args)...));
    }

    void push_front(std::shared_ptr<T> v) {
      std::unique_ptr<Node> new_node(new Node(v));
      std::lock_guard<std::mutex> lock(head_.m);
      new_node->next = std::move(head_.next);
      head_.next = std::move(new_node);
    }

    template<typename Function>
    void for_each(Function f) {
      Node* current = &head_;
      std::unique_lock<std::mutex> lock(head_.m);
      while(Node* const next = current->next.get()) {
        std::unique_lock<std::mutex> next_lock(next->m);
        lock.unlock();
        f(std::ref(*(next->data)));
        current = next;
        lock = std::move(next_lock);
      }
    }

    template<typename Predicate>
    std::shared_ptr<T> pop_first_if(Predicate p) {
      // start searching from the head
      Node* current = &head_;
      std::unique_lock<std::mutex> lock(head_.m);
      while (Node* const next = current->next.get()) {
        std::unique_lock<std::mutex> next_lock(next->m);
        if(p(std::ref(*(next->data)))) {
          // Grab desired data
          std::shared_ptr<T> out = next->data;

          // remove corresponding node
          std::unique_ptr<Node> old_next = std::move(current->next);
          current->next = std::move(old_next->next);
          next_lock.unlock();

          return out;
        }
        else {
          current = next;
          lock.unlock();
          // takes care of unlocking
          lock = std::move(next_lock);
        }
      }
      return std::shared_ptr<T>();
    }

    template <typename Predicate>
    void remove_if(Predicate p) {
      Node* current = &head_;
      std::unique_lock<std::mutex> lock(head_.m);
      while (Node* const next = current->next.get()) {
        std::unique_lock<std::mutex> next_lock(next->m);
        if (p(std::ref(*next->data))) {
          // remove
          std::unique_ptr<Node> old_next = std::move(current->next);
          current->next = std::move(next->next);
          next_lock.unlock();
        }
        else {
          // move on
          current = next;
          // This takes care of unlocking
          lock = std::move(next_lock);
        }
      }
    }
};

}

#endif
