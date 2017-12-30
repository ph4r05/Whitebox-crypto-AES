//
// Created by Dusan Klinec on 30.12.17.
//

#include <gtest/gtest.h>
#include "Commons.h"
#include "../WBAESGenerator.h"

GF2X initDefaultModulus(){
    GF2X gf2x = GF2XFromLong(0x11B, 9);
    GF2E::init(gf2x);
    return gf2x;
}

