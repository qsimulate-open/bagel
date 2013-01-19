//
// BAGEL - Parallel electron correlation program.
// Filename: serverflush.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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

#ifndef __SRC_UTIL_SERVERFLUSH_H
#define __SRC_UTIL_SERVERFLUSH_H

// TODO until GCC fixes this bug
#define _GLIBCXX_USE_NANOSLEEP
#include <atomic>
#include <thread>
#include <src/util/constants.h>

namespace bagel {

// implements a server that periodically flushes
class ServerFlush {
  protected:
    std::atomic<bool> thread_alive_;

    virtual void flush() = 0;

    std::shared_ptr<std::thread> server_;

    void periodic() {
      while (thread_alive_) {
        flush();
        std::this_thread::sleep_for(sleeptime__); 
      }
    }

    std::mutex mut_;

  public:
    ServerFlush() : thread_alive_(false) {
    }

    ~ServerFlush() {
      turn_off();
    }

    void turn_on() {
      std::lock_guard<std::mutex> lock(mut_);
      assert(!thread_alive_);
      thread_alive_ = true;
      server_ = std::shared_ptr<std::thread>(new std::thread(&ServerFlush::periodic, this));
    }

    void turn_off() {
      std::lock_guard<std::mutex> lock(mut_);
      if (thread_alive_) {
        thread_alive_ = false;
        server_->join();
      }
    }
};

}

#endif
