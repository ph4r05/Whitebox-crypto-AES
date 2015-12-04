//
// Created by Dusan Klinec on 04.12.15.
//

#include "base.h"
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>

int phrand(){
    static boost::random_device rd;
    static boost::random::uniform_int_distribution<int> dis;
    return dis(rd);
}
