/*! \file meshlist.h
    \brief アルゴンに対して、分子動力学シミュレーションを行うクラスの宣言

    Copyright ©  2017 @dc1394 All Rights Reserved.
    This software is released under the BSD 2-Clause License.
*/

#ifndef _MESHLIST_H_
#define _MESHLIST_H_

#pragma once

#include "systemparam.h"
#include <cstdint>          // for std::int32_t

namespace moleculardynamics {
    //! A class.
    /*!
        メッシュリストクラス    
    */
    class MeshList final {
    private:
        double mesh_size;
        double periodiclen_;
        std::int32_t m;
        std::int32_t number_of_mesh;
        std::vector<std::int32_t> count;
        std::vector<std::int32_t> indexes;
        std::vector<std::int32_t> sorted_buffer;
        void search(std::int32_t index, SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs);
        void search_other(std::int32_t id, std::int32_t ix, std::int32_t iy, std::int32_t iz, SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs);
    public:
        explicit MeshList(double periodiclen);
        void make_pair(SystemParam::myatomvector & atoms, SystemParam::mypairvector & pairs);
        void set_number_of_atoms(std::size_t pn) { sorted_buffer.resize(pn); }
    };
}

#endif	// _MESHLIST_H_
